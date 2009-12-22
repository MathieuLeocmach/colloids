#include "../graphicParticles.hpp"

using namespace std;
using namespace cimg_library;

/**
    \brief export the xy slices of a 3D image to tiff files.
*/
template <typename T>
void toTiff(const CImg<T> &img, const string &outputPath)
{
    ostringstream oss;
    oss << img.depth;
    const size_t nbZero = string(oss.str()).size();
    cimg_for2Z(img,z)
    {
        ostringstream osz;
        osz << z;
        string zstr(osz.str());
        for(size_t i = zstr.length();i<nbZero;++i)
        {
            zstr = "0"+ zstr;
        }
        img.get_shared_plane(z).save((outputPath+"_z"+zstr+".tif").c_str());
    }
}

int main(int argc, char ** argv)
{
    if(argc<4)
    {
        cout << "Syntax : data2img [path]filename radius outputPath [z_blur [%noise [type of noise]]]" << endl;
        cout << " convert data to a stack of images (~reverse tracking)" << endl;
        return 1;
    }
    const string filename(argv[1]);
    //const string inputPath = filename.substr(0,filename.find_last_of("."));
    float radius;
    sscanf(argv[2],"%f",&radius);
    const string outputPath(argv[3]);
    size_t tNoise;

    GraphicParticles Centers(filename,radius);

    const double mul= 256.0/Centers.getMinDim();
    Centers*=mul;

    CImg<unsigned char> img = Centers.getRepresentation((unsigned char)255,(unsigned char)0);
    if(argc>4)
    {
        float zblur;
        sscanf(argv[4],"%f",&zblur);
        img.blur(0,0,zblur);
    }
    if(argc>5)
    {
        float pNoise;
        sscanf(argv[5],"%f",&pNoise);
        if(argc>6)
        {
            sscanf(argv[6],"%u",&tNoise);
            img.noise(-pNoise,tNoise);
        }
        else img.noise(-pNoise);
    }

    img.display();

    toTiff(img,outputPath);

    return EXIT_SUCCESS;
}
