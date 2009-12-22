#include "../indexedParticles.hpp"

using namespace std;

void addToOrientation(vector< valarray<double> >&orient,const valarray<double>&diff,const double&scale)
{
    for(size_t i=0;i<3;++i)
    {
        orient[(size_t)(orient.size()/2.0+diff[i]*scale/2.0)][i]++;
    }
}


void saveOrient(const vector< valarray<double> > &orient,const string &filename,const float rscale)
{
     cout << "export to " << filename << endl;
     ofstream output(filename.c_str(), ios::out | ios::trunc);
     if(output)
     {
       //header line
       output << "r\tx\ty\tz" << endl;

       for(size_t r=0; r<orient.size();++r)
       {
         output << (r-(float)orient.size()/2.0)/rscale;
         for(size_t i=0; i<3;++i)
           output << "\t" << orient[r][i];
         output << endl;
       }
       output.close();
     }
     else cout << " cannot open the file";
}


int main(int argc, char ** argv)
{
    if(argc<4)
    {
        cout << "Syntax : orientation [path]filename radius NbOfBins" << endl;
        return 1;
    }

	cout << "\t\t\t**************************************" << endl;
    cout << "\t\t\t*** Orientation in the first shell ***" << endl;
    cout << "\t\t\t**************************************" << endl << endl;

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    float radius;
    sscanf(argv[2],"%f",&radius);
    size_t Nbins;
    sscanf(argv[3],"%u",&Nbins);


    /** \brief orientation data container */
    vector< valarray<double> > orient (Nbins, valarray<double>(0.0,3));

    indexedParticles Centers(filename,radius);
    cout << Centers.size() << " particles ... ";
    //The cutoff is at the end of the first peak of the g(r)
    Centers.getHistogram(orient,1.33,addToOrientation);
    cout << " done !" << endl;

    saveOrient(orient,inputPath + "_orientation.txt",Nbins/(2*1.33));
    return EXIT_SUCCESS;
}
