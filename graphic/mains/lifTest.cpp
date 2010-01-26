/**
    Copyright 2008,2009 Mathieu Leocmach

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
*/
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
