/**
    Copyright 2008,2009,2010 Mathieu Leocmach

    This file is part of Colloids.

    Colloids is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Colloids is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.
**/

#include "../dynamicParticles.hpp"
#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
	try
    {
		if(argc<5)
		{
			cerr<<"syntax: linkboo [path]filename token delta_t span [offset=0]"<<endl;
			return EXIT_FAILURE;
		}

		const string filename(argv[1]), token(argv[2]);
		const double delta_t = atof(argv[3]);
		const size_t span = atoi(argv[4]);
		const size_t offset = (argc>5)?atoi(argv[5]):0;

		cout<<filename.substr(filename.find_last_of("/\\")+1)<<endl;

		//create the needed file series
		FileSerie datSerie(filename, token, span, offset),
				bondSerie = datSerie.changeExt(".bonds"),
				qlmSerie = datSerie.changeExt(".qlm"),
				cloudSerie = datSerie.changeExt(".cloud"),
				cgCloudSerie = datSerie.addPostfix("_space", ".cloud");

		cout<<"load ..."<<endl;
		//load all files in memory with default radius of 1.0
		//by the way, check file existence
		boost::ptr_vector<Particles> positions(span);
		for(size_t t=0; t<span; ++t)
			positions.push_back(new Particles(datSerie%t));

		//spatially index each frame
		cout<<"index ..."<<endl;
		for_each(positions.begin(), positions.end(), mem_fun_ref(&Particles::makeRTreeIndex));

		//get averaged g(r)
		cout<<"g(r) in "<<(datSerie.head()+".rdf")<<endl;
		vector<double> total_g(200, 0.0), g;
		boost::progress_display show_pr(span);
		for(size_t t=0; t<span; ++t)
		{
			g = positions[t].getRdf(200,15.0);
			transform(
				g.begin(), g.end(),
				total_g.begin(), total_g.begin(),
				plus<double>()
				);
			++show_pr;
		}
		transform(
			total_g.begin(), total_g.end(), total_g.begin(),
			bind2nd(divides<double>(), span)
			);
		ofstream rdfFile((datSerie.head()+".rdf").c_str(), ios::out | ios::trunc);
		rdfFile << "#r\tg(r)"<<endl;
		for(size_t r=0; r<total_g.size(); ++r)
			rdfFile<< r/200.0*15.0 <<"\t"<< total_g[r]<<"\n";

		//get bondlength and radius from the first minimum of g(r)
		//the loop is here only to get rid of possible multiple centers at small r
		vector<double>::iterator first_peak = total_g.begin();
		size_t first_min;
		do
		{
			first_peak = max_element(total_g.begin(),total_g.end());
			first_min = distance(total_g.begin(), min_element(first_peak,total_g.end()));
			//cout<<"first_peak="<<distance(total_g.begin(), first_peak)<<" first_min="<<first_min<<" ... ";
		}
		while(total_g[first_min]==0.0);
		const double bondLength = first_min/200.0*15.0, radius = bondLength / 1.3;
		cout<<"radius="<<radius<<endl;

		//treat each file
		boost::progress_display show_progress(span);
		for(size_t t=0; t<span; ++t)
		{
			//create neighbour list and export bonds
			positions[t].makeNgbList(bondLength);
			BondList bonds = positions[t].getBonds();
			ofstream bondFile((bondSerie%t).c_str(), ios::out | ios::trunc);
			for(deque<pair<size_t, size_t> >::const_iterator b=bonds.begin(); b!= bonds.end();++b)
				bondFile<<b->first<<" "<<b->second<<"\n";
			bondFile.close();

			//select the particles further than the bond length from the boundaries
			set<size_t> inside = positions[t].selectInside(bondLength),
                secondInside = positions[t].selectInside(2.0*bondLength);

			//calculate and export qlm
			vector<BooData> qlm, qlm_cg;
			positions[t].getBOOs(inside, qlm);
			positions[t].getCgBOOs(secondInside, qlm, qlm_cg);
			ofstream qlmFile((qlmSerie%t).c_str(), ios::out | ios::trunc);
			copy(
				qlm_cg.begin(), qlm_cg.end(),
				ostream_iterator<BooData>(qlmFile,"\n")
				);

			//calculate and export invarients
			ofstream cloudFile((cloudSerie%t).c_str(), ios::out | ios::trunc);
			cloudFile<<"#Q4\tQ6\tW4\tW6"<<endl;
			transform(
				qlm.begin(), qlm.end(),
				ostream_iterator<string>(cloudFile,"\n"),
				cloud_exporter()
				);
			cloudFile.close();

			ofstream cloud_cgFile((cgCloudSerie%t).c_str(), ios::out | ios::trunc);
			cloud_cgFile<<"#Q4\tQ6\tW4\tW6"<<endl;
			transform(
				qlm_cg.begin(), qlm_cg.end(),
				ostream_iterator<string>(cloud_cgFile,"\n"),
				cloud_exporter()
				);
			cloud_cgFile.close();

			//update radius
			positions[t].radius = radius;
			++show_progress;
		}

		//link and save trajectories
		DynamicParticles(positions, radius, delta_t).save(
			datSerie.head()+".traj", filename.substr(filename.find_last_of("/\\")+1), token, offset, span
			);
	}
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

