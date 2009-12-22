#include "pv.hpp"

using namespace std;

/** \brief input to a binary file */
ostream& operator<< (ostream& out, const atom &ato )
{
    out.write((char*)&ato.attr,sizeof(unsigned char));
    out.write((char*)&ato.x,sizeof(long));
    out.write((char*)&ato.y,sizeof(long));
    out.write((char*)&ato.z,sizeof(long));
    return out;
}
/** \brief import from a binary file */
istream& operator>> (istream& is, atom& ato )
{
    is.read((char*)&ato.attr,sizeof(unsigned char));
    is.read((char*)&ato.x,sizeof(long));
    is.read((char*)&ato.y,sizeof(long));
    is.read((char*)&ato.z,sizeof(long));
    return is;
}


/** \brief constructor from a pv binary file */
pv::pv(const string &filename) : deque<Particles>(0,Particles(0,0.0))
{
    ifstream file(filename.c_str(), ios::in | ios::binary);
    if(!file)
		throw invalid_argument("No such file as "+filename);

	char ident[32];
	file.read((char*)&ident,32*sizeof(char));
	name=string(ident);
	size_t ndata =0;
	if(name.compare(0,32,"PV-01 /Shoji-Maruyama Laboratory")==0)
	{
		header16 h;
		file.read((char*)&h,sizeof(header16));
		moltype=h.moltype;
		nmol=h.nmol;
		nbond=h.nbond;
		//h.ndata=(unsigned short)size();
		bb.edges[0].first = h.vlx1;
		bb.edges[0].second= h.vlx2;
		bb.edges[1].first = h.vly1;
		bb.edges[1].second= h.vly2;
		bb.edges[2].first = h.vlz1;
		bb.edges[2].second= h.vlz2;
		time0=h.time0;
		dt=h.dt;
		ndata=h.ndata;
	}
	else
	{
		if(name.compare(0,32,"PV-32 /Shoji-Maruyama Laboratory")==0)
		{
			header32 h;
			file.read((char*)&h,sizeof(header32));
			moltype=h.moltype;
			nmol=h.nmol;
			nbond=h.nbond;
			//h.ndata=(unsigned short)size();
			bb.edges[0].first = h.vlx1;
			bb.edges[0].second= h.vlx2;
			bb.edges[1].first = h.vly1;
			bb.edges[1].second= h.vly2;
			bb.edges[2].first = h.vlz1;
			bb.edges[2].second= h.vlz2;
			time0=h.time0;
			dt=h.dt;
			ndata=h.ndata;
		}
		else
		{
			cerr << filename <<" is not a valid PV file" << endl;
			return;
		}
	}

	file.read((char*)&diam,sizeof(double));

	atom ato;
	valarray<double> p(0.0,3);
	labels.assign(ndata,std::map<size_t,unsigned char>());
	for(size_t i=0;i<ndata;++i)
	{
		Particles ps(0,0.0);
		ps.bb=bb;
		for(int j=0; j<nmol;++j)
		{
			file>>ato;
			if(ato.x!=-10000 || ato.y!=-10000 || ato.z!=-10000)
			{
				p[0]=(double)ato.x/100.0;
				p[1]=(double)ato.y/100.0;
				p[2]=(double)ato.z/100.0;
				ps.push_back(p);
				labels[i].insert(labels[i].end(),make_pair(j,ato.attr));
				//ps.labels.push_back(ato.attr);
			}
		}
		push_back(ps);
	}
	file.close();

    return;
}


/** \brief appends a frame at the end */
void pv::push_back(const Particles &ps,const std::map<size_t,unsigned char> *lab)
{
    deque<Particles>::push_back(ps);
    if(lab)
		labels.resize(size(),*lab);
	else
		labels.resize(size(),std::map<size_t,unsigned char>());
    bb.stretch(ps.bb);
    //dt++;
    nmol = max(nmol,(long)ps.size());
}
/** \brief concatenates with another animation */
void pv::operator<<(const pv &a)
{
    const size_t old_end = size();
    for(const_iterator ps=a.begin();ps!=a.end();++ps)
        push_back(*ps);
	labels.resize(size());
	for(size_t i=0;i<a.labels.size();++i)
		labels[old_end+i] = a.labels[i];
}

/** \brief exports to binary PV files */
void pv::exportToPV(const string &filename)
{
    cout << "export to " << filename << endl;

    ofstream fileoutput(filename.c_str(), ios::out | ios::trunc | ios::binary);
    if(fileoutput)
    {
        //PV header
        char ident[32];
        name.copy(ident,32);
        fileoutput.write((char*)&ident,32*sizeof(char));
        header32 h;
        h.moltype=moltype;
        h.nmol=nmol;
        h.nbond=nbond;
        h.ndata=(long)size();
        h.vlx1=bb.edges[0].first;
        h.vlx2=bb.edges[0].second;
        h.vly1=bb.edges[1].first;
        h.vly2=bb.edges[1].second;
        h.vlz1=bb.edges[2].first;
        h.vlz2=bb.edges[2].second;
        h.time0=time0;
        h.dt=dt;

        fileoutput.write((char*)&h,sizeof(header32));
        fileoutput.write((char*)&diam,sizeof(double));

        atom atomNull,ato;
        atomNull.attr=0;
        atomNull.x = (short)-10000;
        atomNull.y = (short)-10000;
        atomNull.z = (short)-10000;
        ato.attr=1;
        //struct atom * atombuf=(struct atom *)malloc(nmol*sizeof(struct atom));

        map<size_t,unsigned char>::const_iterator lab;
        for(size_t t=0;t<size();++t)
        {
            for(size_t p=0;p<at(t).size();++p)
            {
                ato.x = (short)(at(t)[p][0]*100);
                ato.y = (short)(at(t)[p][1]*100);
                ato.z = (short)(at(t)[p][2]*100);

                lab = labels[t].find(p);
				if(lab!=labels[t].end())
					ato.attr=(*lab).second;
                else
					ato.attr=1;
                fileoutput << ato;
            }
            for(size_t p=at(t).size();p<(size_t)nmol;++p)
            {
                fileoutput <<atomNull;
            }

            //fileoutput.write(output.c_str(),nmol*sizeof(atom));
            //fwrite(atombuf,sizeof(struct atom),nmol,fileoutput.rdbuf());
        }
        fileoutput.close();
    }
    else cout << " cannot open the file";
}
