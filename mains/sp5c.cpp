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
		if(argc<2)
		{
			cerr<<"syntax: sp5c [path]filename"<<endl;
			return EXIT_FAILURE;
		}

		const string filename(argv[1]),
			ext = filename.substr(filename.find_last_of(".")),
			path = filename.substr(0, filename.find_last_of("/\\")+1),
			inputPath = filename.substr(0,filename.find_last_of("."));

        //load data
        Particles parts(filename);
        BondSet bonds = loadBonds(inputPath+".bonds");
        parts.makeNgbList(bonds);
        //get the pair of particles having 5 common neighbours
        list<size_t> sp5c;
        //count the number of 1551 paris a given particle is part of.
        //The central particle of a perfect icosahedron is part of 12 1551 pairs
        vector<size_t> spindle(parts.size(),0);
        for(BondSet::const_iterator b=bonds.begin(); b!=bonds.end(); ++b)
        {
            //stringstream s;
            //find the common neighbours of the two extremities of the bond
            list<size_t> common;
            set_intersection(
                parts.getNgbList()[b->low()].begin(), parts.getNgbList()[b->low()].end(),
                parts.getNgbList()[b->high()].begin(), parts.getNgbList()[b->high()].end(),
                back_inserter(common)
                );
            if(common.size()==5)
            {
                map<size_t, list<size_t> > ringngb;
                bool is_ring = true;
                for(list<size_t>::const_iterator c=common.begin(); c!=common.end(); ++c)
                {
                    set_intersection(
                        parts.getNgbList()[*c].begin(), parts.getNgbList()[*c].end(),
                        common.begin(), common.end(),
                        back_inserter(ringngb[*c])
                        );
                    is_ring &= (ringngb[*c].size()==2);
                    //s<<*c<<" linked to ";
                    //copy(ringngb[*c].begin(), ringngb[*c].end(), ostream_iterator<size_t>(s, "\t"));
                    //s<<"\n";
                }
                if(!is_ring) continue;


                //p=0;
                //ringngb[0].front()

                sp5c.push_back(b->low());
                sp5c.push_back(b->high());
                spindle[b->low()]++;
                spindle[b->high()]++;
                //sort the common neighbours in order to have a convex pentagon
                sp5c.push_back(common.front());
                size_t last = ringngb[sp5c.back()].back(), previous=last;
                do
                {
                    //s<<previous<<"->"<<sp5c.back()<<"->";
                    if(ringngb[sp5c.back()].front()!=previous)
                    {
                        previous = sp5c.back();
                        sp5c.push_back(ringngb[previous].front());
                        //s<<sp5c.back()<<"\n";
                    }
                    else
                    {
                        previous = sp5c.back();
                        sp5c.push_back(ringngb[previous].back());
                        //s<<sp5c.back()<<"\n";
                    }
                }while(previous!=last);
                sp5c.pop_back();
                //cout<<s.str()<<endl;
                //exit(0);

                //sp5c.splice(sp5c.end(), common);
            }

        }
        assert(sp5c.size()%7==0);
        cout<<sp5c.size()/7<<" pairs with 5 common neighbours in a ring"<<endl;

        cout<<"bond orientational order"<<endl;
        vector<BooData> qlm(parts.size());
        parts.getBOOs(qlm);

        cout<<"invariants"<<endl;
        vector< boost::array<double, 8> > invariants(qlm.size());
        for(size_t p=0; p<qlm.size();++p)
            for(size_t l=2; l<=10; l+=2)
                qlm[p].getInvarients(l, invariants[p][l/2-1], invariants[p][l/2+2]);

        //export to vtk
        cout<<"export vtk"<<endl;

        ofstream out((inputPath+"_ico.vtk").c_str());

        parts.toVTKstream(out, "Icosahedral order");

        toVTKstream(out, bonds);

        out<<"POINT_DATA "<<invariants.size()<<"\n";
        for(int i=0;i<8;++i)
        {
            out<<"SCALARS "<<((i/4)?"W":"Q")<<((i%4)*2+4)<<" float 1\n";
            out<<"LOOKUP_TABLE default\n";
            for(size_t p=0; p<invariants.size(); ++p)
                out<<invariants[p][i]<<"\n";
        }
        out<<"SCALARS spindle int 1\n";
        out<<"LOOKUP_TABLE default\n";
        copy(spindle.begin(), spindle.end(), ostream_iterator<size_t>(out, " "));


        out.close();

    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
