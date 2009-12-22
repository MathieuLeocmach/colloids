#include "../lifFile.hpp"
#include "CImg.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "Syntax : liftest [path]filename" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    LifFile lif(filename);
    cout << "LIF version "<<lif.LifVersion << endl;
    size_t serie = lif.chooseSerie(),frame=0;
    if(lif.Dimensions[serie]->size()>3)
        do
        {
            cout<<"Chose a frame between 0 and "<<lif.Dimensions[serie]->at(3)->NumberOfElements-1<<": ";
            cin>>frame;
        }
        while(lif.Dimensions[serie]->at(3)->NumberOfElements<frame);
    CImg<unsigned char>img;
    switch (lif.Dimensions[serie]->size())
    {
        case 0:
            cerr << "Serie "<<serie<<" doesn't exist or is empty"<<endl;
            return EXIT_FAILURE;
        case 1:
            img.assign(
                lif.Dimensions[serie]->at(0)->NumberOfElements);
            break;
        case 2:
            img.assign(
                lif.Dimensions[serie]->at(0)->NumberOfElements,
                lif.Dimensions[serie]->at(1)->NumberOfElements);
            break;
        default :
            img.assign(
                lif.Dimensions[serie]->at(0)->NumberOfElements,
                lif.Dimensions[serie]->at(1)->NumberOfElements,
                lif.Dimensions[serie]->at(2)->NumberOfElements);
    }
    cout<<"image constructed"<<endl;
    lif.fill3DBuffer(img.ptr(),serie,frame);
    img.display();

    return EXIT_SUCCESS;
}
