/**
 * \file tracker.cpp
 * \brief Implement classes to track particles in 3D image series
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 7 april 2009
 *
 * Object oriented implementation of the code elaborated in 2008.
 * Old auto-configuration routine are commented at the end but not implemented this time
 * Based on Crocker & Grier algorithm
 *
 */


#include "tracker.hpp"
#include <boost/progress.hpp>

using namespace std;
using namespace cimg_library;

/** \brief constructor from a serie of a lif file */
LifTracker::LifTracker(std::auto_ptr_ref<LifFile> l,const size_t &s,const size_t &ch,const unsigned flags)
{
    lif = l;
    serie = s;
    if(ch>0 && lif->Channels[serie]->size()<=ch)
    {
        ostringstream os;
        os<<lif->Names[serie]<<" has only "<<lif->Channels[serie]->size()<<" channels";
        throw invalid_argument(os.str());
    }
    channel = ch;
    switch (lif->Dimensions[serie]->size())
    {
        case 0:
            {
                ostringstream os;
                os << "Serie "<<serie<<" doesn't exist or is empty";
                throw invalid_argument(os.str());
            }
            break;
        case 1:
            image.assign(
                lif->Dimensions[serie]->at(0)->NumberOfElements);
            data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)
                * lif->Dimensions[serie]->at(0)->NumberOfElements);
            forward_plan = fftwf_plan_dft_1d(lif->Dimensions[serie]->at(0)->NumberOfElements, data, data, FFTW_FORWARD, flags);
            backward_plan = fftwf_plan_dft_1d(lif->Dimensions[serie]->at(0)->NumberOfElements, data, data, FFTW_BACKWARD, flags);
            break;
        case 2:
            image.assign(
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements);
            data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)
                * lif->Dimensions[serie]->at(0)->NumberOfElements
                * lif->Dimensions[serie]->at(1)->NumberOfElements);
            forward_plan = fftwf_plan_dft_2d(
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements,
                data, data, FFTW_FORWARD, flags);
            backward_plan = fftwf_plan_dft_2d(
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements,
                data, data, FFTW_BACKWARD, flags);
            break;
        default :
            const size_t nbPix = lif->Dimensions[serie]->at(0)->NumberOfElements
                * lif->Dimensions[serie]->at(1)->NumberOfElements
                * lif->Dimensions[serie]->at(2)->NumberOfElements;

            data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nbPix);
            cout<<"fftw data allocated"<<endl;
            forward_plan = fftwf_plan_dft_3d(
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements,
                lif->Dimensions[serie]->at(2)->NumberOfElements,
                data, data, FFTW_FORWARD, flags);
            cout<<"forward plan ok"<<endl;
            backward_plan = fftwf_plan_dft_3d(
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements,
                lif->Dimensions[serie]->at(2)->NumberOfElements,
                data, data, FFTW_BACKWARD, flags);
            cout<<"backward plan ok"<<endl;
            //imgdata = (float*)fftwf_malloc(sizeof(float) * nbPix);
            image.assign(
                //imgdata,
                lif->Dimensions[serie]->at(0)->NumberOfElements,
                lif->Dimensions[serie]->at(1)->NumberOfElements,
                lif->Dimensions[serie]->at(2)->NumberOfElements);//,
                //1,true);
            cout<<"image allocated"<<endl;
    }
    //dilated.assign(((float*)data)+image.size(),image.dimx(),image.dimy(),image.dimz(),1,true);
    cout<<"dilate allocated"<<endl;
    Zratio = lif->getVoxelSize(serie,2)/lif->getVoxelSize(serie,0);
}

/** @brief ~LifTracker destructor  */
 LifTracker::~LifTracker()
{
    fftwf_destroy_plan(forward_plan);
	fftwf_destroy_plan(backward_plan);
	fftwf_free(data);
}

/** \brief load a 3D stack from the current lif serie*/
void LifTracker::load3D(const size_t &time)
{
    //CImg has non interleaved channels wheras lif file data's channels are interleaved.
    //If their is more than one channel we need to use an intermediate image and permute its axes
    if(lif->Channels[serie]->size()>1)
    {
        CImg<unsigned char> colorful(lif->Channels[serie]->size(),image.width,image.height,image.depth);
        lif->fill3DBuffer(colorful.ptr(),serie,time);
        image=colorful.permute_axes("yzvx").channel(channel);
    }
    else
    {
        CImg<unsigned char> temp(image,"xyzv");
        lif->fill3DBuffer(temp.ptr(),serie,time);
        image = temp;
    }
    if(!quiet)  cout<<"t="<<time<<" image loaded"<<endl;
}


/** \brief constructor from a config file (to read tiff) */
FileSerieTracker::FileSerieTracker(TokenTree &tt,const double &Zra,const size_t &xsize, const size_t &ysize,const size_t &zsize,const size_t &ch,unsigned flags)
{
    fileSerie = &tt;
    image.assign(xsize,ysize,zsize);
    //dilated.assign(image,"xyzv");
    data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * xsize * ysize * zsize);
    forward_plan = fftwf_plan_dft_3d(xsize,ysize,zsize, data, data, FFTW_FORWARD, flags);
    backward_plan = fftwf_plan_dft_3d(xsize,ysize,zsize, data, data, FFTW_BACKWARD, flags);
    Zratio = Zra;
    channel = ch;
}

/** \brief load a 3D stack from file serie */
void FileSerieTracker::load3D(const size_t &time)
{
    if(image.dimz()==1)
    {
        if(fileSerie->nbdigit==0)
            image=CImg<unsigned char>(fileSerie->getPattern().c_str()).channel(channel);
        else
            image=CImg<unsigned char>(((*fileSerie) %time).str().c_str()).channel(channel);
    }
    else
    {

        if(!quiet)
        {
            cout <<"t="<<time<<" loading image files"<<endl;
            boost::progress_display show_progress(image.dimz());
            for(size_t z=0;z<(size_t)image.dimz();++z)
            {
                ++show_progress;
                image.get_shared_plane(z)=CImg<unsigned char>(((*fileSerie) %time % z).str().c_str()).channel(channel);
            }
        }
        else
            for(size_t z=0;z<(size_t)image.dimz();++z)
                image.get_shared_plane(z)=CImg<unsigned char>(((*fileSerie) %time % z).str().c_str()).channel(channel);
    }
}
/** @brief make the band pass filter mask in Fourier space
    \param radiusMin Real space mimimum radius in (x,y) plane
    \param radiusMax Real space maximum radius in (x,y) plane
    \param zradiusMin Real space mimimum radius in z
    \param zradiusMax Real space maximum radius in z
*/
void Tracker::makeBandPassMask(const double &radiusMin, const double &radiusMax, const double &zradiusMin, const double &zradiusMax)
{
    if(!quiet) cout << "making the band-pass filter ... ";
    const double freqMax = image.dimx()/2.0/radiusMin,
              freqMin =image.dimx()/2.0/radiusMax,
              freqzMax = image.dimz()/2.0/zradiusMin,
              freqzMin =image.dimz()/2.0/zradiusMax;

    FFTmask.assign(image,"xyzv",0);
    drawEllipsoid(FFTmask,image.dimx()/2,image.dimy()/2,image.dimz()/2,freqMax,freqzMax,true);
    drawEllipsoid(FFTmask,image.dimx()/2,image.dimy()/2,image.dimz()/2,freqMin,freqzMin,false);

    //The spectrum output of FFT is centered on (0,0), not in the center of the image
    //It is lighter to (translate mask+apply mask) rather than (translate spectrum+apply mask+translate back spectrum)
    FFTmask.translate(image.dimx()/2,image.dimy()/2,image.dimz()/2,0,2);
}


/**
    \brief Fourier transform, apply Mask, inverse FFT and normalize the result between 0 and 255.
    The heaviest part of the algorithm
*/
void Tracker::FFTapplyMask()
{
    if(!quiet) cout << "fill in data ... ";
    //fill in the content of image into the fftw data (row major)
    const unsigned int w = image.width, wh = w*image.height, whd = wh*image.depth;
    float *ptrr = image.ptr(0,0,0,0), *ptrd = (float*)data;
    cimg_forX(image,x)
    {
        cimg_forY(image,y)
        {
            cimg_forZ(image,z)
            {
                *(ptrd++) = *ptrr;
                *(ptrd++) = 0.0;
                ptrr+=wh;
            }
            ptrr-=whd-w;
        }
        ptrr-=wh-1;
    }

    if(!quiet) cout << "FFT ... ";
    fftwf_execute(forward_plan);

    if(!quiet) cout << "applying filter ... ";
    bool *ptrm = FFTmask.ptr(0,0,0,0);
    ptrd = (float*)data;
    cimg_forX(FFTmask,x)
    {
        cimg_forY(FFTmask,y)
        {
            cimg_forZ(FFTmask,z)
            {
                *(ptrd++) *= *ptrm;
                *(ptrd++) *= *ptrm;
                ptrm+=wh;
            }
            ptrm-=whd-w;
        }
        ptrm-=wh-1;
    }

    if(!quiet) cout << "inverse FFT ... ";
    fftwf_execute(backward_plan);

    if(!quiet) cout << "get data back ... ";
    ptrd = (float*)data;
    ptrr = image.ptr(0,0,0,0);
    cimg_forX(image,x)
    {
        cimg_forY(image,y)
        {
            cimg_forZ(image,z)
            {
                *ptrr = *(ptrd++);
                ptrd++;
                ptrr+=wh;
            }
            ptrr-=whd-w;
        }
        ptrr-=wh-1;
    }
    if(!quiet) cout << "normalize ... ";
    image.normalize(0,255);
}

/** \brief  Dilate the image by a structuring element.
    It is a copy of CImg code but with no new memory allocation
*/
CImg<float>& get_dilate(const CImg<bool>& mask, const CImg<float>&src,CImg<float>&dest) {
    if (src.is_empty()) return dest;
    if (!mask || mask.dim!=1)
        throw CImgArgumentException("CImg<%s>::dilate() : Specified mask (%u,%u,%u,%u,%p) is not a scalar image.",
                                src.pixel_type(),mask.width,mask.height,mask.depth,mask.dim,mask.data);
    const int
        mx2 = mask.dimx()/2, my2 = mask.dimy()/2, mz2 = mask.dimz()/2,
        mx1 = mx2 - 1 + (mask.dimx()%2), my1 = my2 - 1 + (mask.dimy()%2), mz1 = mz2 - 1 + (mask.dimz()%2),
        mxe = src.dimx() - mx2, mye = src.dimy() - my2, mze = src.dimz() - mz2;
    const int v=0;
    for (int z = mz1; z<mze; ++z) for (int y = my1; y<mye; ++y) for (int x = mx1; x<mxe; ++x)
    {
        float max_val = cimg::type<float>::min();
        for (int zm = -mz1; zm<=mz2; ++zm) for (int ym = -my1; ym<=my2; ++ym) for (int xm = -mx1; xm<=mx2; ++xm)
        {
            const float cval = src(x+xm,y+ym,z+zm,v);
            if (mask(mx1+xm,my1+ym,mz1+zm) && cval>max_val) max_val = cval;
        }
        dest(x,y,z,v) = max_val;
    }
    cimg_forYZV(src,y,z,v)
        for (int x = 0; x<src.dimx(); (y<my1 || y>=mye || z<mz1 || z>=mze)?++x:((x<mx1-1 || x>=mxe)?++x:(x=mxe)))
        {
            float max_val = cimg::type<float>::min();
            for (int zm = -mz1; zm<=mz2; ++zm) for (int ym = -my1; ym<=my2; ++ym) for (int xm = -mx1; xm<=mx2; ++xm)
            {
                const float cval = src._atXYZ(x+xm,y+ym,z+zm,v);
                if (mask(mx1+xm,my1+ym,mz1+zm) && cval>max_val) max_val = cval;
            }
            dest(x,y,z,v) = max_val;
        }
  return dest;
}

/**
    \brief The main tracking function.
*/
GraphicParticles Tracker::trackXYZ(const double &threshold)
{
    if(view) image.display("Original stack");

    if(!quiet) cout << "band-passing ... ";
    FFTapplyMask();
    if(!quiet) cout << "done!" << endl;
    if(view) image.display("Band Passed");

	//the memory allocated to the FFTW data can be used for other purposes
	CImg<float> dilated(((float*)data)+image.size(),image.dimx(),image.dimy(),image.dimz(),1,true);
	CImg<float> dilated2(((float*)data),image.dimx(),image.dimy(),image.dimz(),1,true);
    //local maxima are found by dilatation of the image (see Peter Lu's PhD thesis)
    //The structuring element is a cube (square), so the dilatation can be separated along each axis
    CImg<bool> dilateMask(3,1,1,1,1);
    if(!quiet) cout << "dilating ... ";
    get_dilate(dilateMask,image,dilated);
    dilateMask.permute_axes("yxzv");
    get_dilate(dilateMask,dilated,dilated2);
    if (image.dimz()>2)
	{
		dilateMask.permute_axes("xzyv");
		get_dilate(dilateMask,dilated2,dilated);
	}
	else
		dilated2.transfer_to(dilated);
	if(view) dilated.display("Dilated");

    /*CImg<bool> dilateMask;
    if (image.dimz()>2)
        dilateMask.assign(3,3,3,1,1);
    else
        dilateMask.assign(3,3,1,1,1);
    if(!quiet) cout << "dilating ... ";
    get_dilate(dilateMask,image,dilated);*/
    if(!quiet) cout << "substracting ... ";
    dilated-=image;
    dilated*=-1;
    if(!quiet) cout << "Find Centers ... ";
    CImg<bool> PixelCenters = dilated.exp().threshold(0.9999);

    // Discard center with intensity < threshold on the band passed image
    //need to use home made implementation to prevent new memory allocation
    if(threshold!=0)
    {
        if(!quiet) cout << "Thresholding ... ";
        for(int i=0;i<image.dimx();++i)
            for(int j=0;j<image.dimy();++j)
                for(int k=0;k<image.dimz();++k)
                    if(image(i,j,k)<threshold)
                        PixelCenters(i,j,k)=false;
    }

    if(!quiet) cout << "done!" << endl;
    if(view) (PixelCenters).display("Found centers");

    if(!quiet) cout << "Vectorize Centers ... ";
    //get the centers in vector format out of the bright pixels
    GraphicParticles Centers(PixelCenters,displayRadius);
    if(!quiet) cout << Centers.size() << " centers found ... done!" << endl;

    if(!quiet) cout << "Sort by decreasing intensity ... ";
    Centers.sortDecreasingIntensity(image);
    if(!quiet) cout << "the smallest intensity is " << Centers.get_neigbourhood(Centers.size()-1,image,1.0).sum() << " ... done!" << endl;

    if(!quiet) cout << "Sub pixel resolution ... ";
    Centers.centroid(image);

    if(!quiet) cout << "Clean Centers' list ... ";
    Centers.clean();
    if(!quiet) cout << Centers.size() << " centers left ... done!" << endl;

    if(view) Centers.displaySpheres(image,150.0f);

    if(image.dimz()>1 && Zratio!=1.0)
    {
        if(!quiet) cout << "Z scaling ... ";
        valarray<double> redim(1.0,3);
        redim[2] = Zratio;
        Centers*=redim;
    }

    if(!quiet) cout << "done!" << endl << endl;

    if(image.dimz()>1)
    {
        if(!quiet) cout << "Calculated volume Fraction : " << Centers.getVF() << endl;
    }
    else
        if(!quiet) cout << "Calculated area Fraction : " << Centers.size() * pow(displayRadius,2.0)*M_PI/(image.dimx()-2.0*displayRadius)/(image.dimy()-2.0*displayRadius) << endl;

    return Centers;
}
/** @brief granulometry
  *
  * @todo: document this function
  */
vector<GraphicParticles*> Tracker::granulometry(const double &radiusMin, const double &radiusMax)
{
	vector<GraphicParticles*> parts((size_t)(radiusMax-radiusMin));
	CImg<bool> PixelCenters(image,"xyz",false);
	CImg<bool> dilateMask;
    if (image.dimz()>2)
        dilateMask.assign(3,3,3,1,1);
    else
        dilateMask.assign(3,3,1,1,1);
	for(size_t s=0;s<parts.size();++s)
	{
		image.blur(radiusMin+s*(radiusMax-radiusMin)/parts.size());
		PixelCenters = (image-image.get_dilate(dilateMask)).exp().threshold(0.9999);
		parts[s] = new GraphicParticles(PixelCenters,radiusMax);
	}
	return parts;
}


/** @brief save FFTW wisdom (plans for fastest FFT)  */
void Tracker::saveWisdom(const std::string & filename)
{
	//no overwritting of the content of an existing file
	loadWisdom(filename);
	FILE * wis = fopen(filename.c_str(),"w+");
	if(wis == NULL)
		throw std::invalid_argument("No such file as "+filename);
	fftwf_export_wisdom_to_file(wis);
	fclose(wis);
}

/** @brief load FFTW wisdom from file
  * \return	false if the file isn't found
  */
bool Tracker::loadWisdom(const std::string & filename)
{
	FILE * wis = fopen(filename.c_str(),"r");
	if(wis == NULL) return false;
	fftwf_import_wisdom_from_file(wis);
	fclose(wis);
	return true;
}


    //}@

    //--------------------------------------
    //
    //! \name Auto-Configuration
    //@{
    //--------------------------------------

/**
    \brief Track centers for a certain value of a given parameter
    \param img Preprocessed image to be tracked
    \param config ConfigFile containing all the tracking parameters
    \param name Name of the parameter to change
    \param value Value of the parameter
    \return Number of tracked particles
*/
/*template<typename T>
int getNbFromParam(const CImg<T> &img,ConfigFile config,const string name, const double value)
{
    config.add(name,value);
    GraphicParticles Centers = trackXYZ(img,config);
    double stepValue;
    config.readInto(stepValue,"stepValue",0.1);
    string outputPath;
    config.readInto<string>(outputPath,"outputPath","");
    ostringstream oss;
    oss << "_t" << (int)(value/stepValue);
    Centers.exportToFile(outputPath+name+oss.str());
    return Centers.size();
}*/


/**
    \brief auto-calibration of the parameters
*/

/*void autoconfig(ConfigFile &config)
{
    //clock_t initial_time = clock ();

    cout << "**************************" << endl;
    cout << "*** AUTO CONFIGURATION ***" << endl;
    cout << "**************************" << endl << endl;

    //read and parse the image name from the configuration file
    string Imagefile,file_prefix,file_sufix,outputPath, time_prefix;
    int nbdigit_z,nbdigit_t;
    if(!config.readInto<string>( Imagefile, "Imagefile" ))
      cerr << "the \"Imagefile\" entry is missing in config.ini";
    parseImageFileName(Imagefile,file_prefix,file_sufix,nbdigit_z,nbdigit_t);

    config.readInto<string>(outputPath,"outputPath","");

    //get the configuration data needed to open files
    int xsize, ysize,zsize,zoffset,tsize,toffset,channel;
    config.readInto(xsize,"xsize");
    config.readInto(ysize,"ysize");
    config.readInto(zsize,"zsize",1);
    config.readInto(zoffset,"zoffset",0);
    config.readInto(tsize,"tsize",1);
    config.readInto(toffset,"toffset",0);
    config.readInto(channel,"channel",0);

    double VoxelWidth,VoxelDepth,Zratio=1;
    if(zsize!=1 && (!config.readInto(VoxelWidth,"VoxelWidth") || !config.readInto(VoxelDepth,"VoxelDepth")) )
      cerr << "the \"VoxelWidth\" or \"VoxelDepth\" entry is missing in config.ini\nLook into the data from the microscope or XML export.";
    else
      Zratio = VoxelWidth / VoxelDepth;


    //Creates a 3D image of the size of each frame of the dataset
    CImg<float> image(xsize,ysize,zsize);

    if(nbdigit_t!=0)
    {
        ostringstream oss;
        oss << toffset;
        string tstr =oss.str();

        for(int j = tstr.length();j<nbdigit_t;j++)
        {
         tstr = "0"+ tstr;
        }

        time_prefix = file_prefix + tstr;
        if(nbdigit_z!=0) time_prefix += "_z";
        cout << endl;
        cout << "\t\t\t*******************" << endl;
        cout << "\t\t\t***time " << tstr << " ***" << endl;
        cout << "\t\t\t*******************" << endl << endl;

        openRoutineZ(image,time_prefix,file_sufix,nbdigit_z,zoffset,channel);
    }
    else
        openRoutineZ(image,file_prefix,file_sufix,nbdigit_z,zoffset,channel);

    //extract from the config file the name of the parameter to be varied
    string nameParam;
    if(!config.readInto(nameParam,"nameParam"))
    {
      cout << "To use the autoconfiguration, you need to fill in a \"nameParam\" parameter into confing.ini" << endl;
      cout << "It indicates the name of the parameter to be auto-configured. For exemple \"nameParam = radius\"." << endl;
    }

    //extract from the config file the bound of the variation
    double minValue,maxValue,stepValue;
    if(!config.readInto(minValue,"minValue"))
      cout << "To use the autoconfiguration, you need to fill in a \"minValue\" parameter into confing.ini" << endl;
    if(!config.readInto(maxValue,"maxValue"))
      cout << "To use the autoconfiguration, you need to fill in a \"maxValue\" parameter into confing.ini" << endl;


    //map of the number of detected particles function of the varying parameter
    map<double,int> Map;

    cout << "\t\ttesting various " << nameParam << " between " << minValue << " and " << maxValue << endl;
    if(!config.readInto(stepValue,"stepValue",0.0))
    {
        Map[minValue]=getNbFromParam(image,config,nameParam,minValue);
        Map[maxValue]=getNbFromParam(image,config,nameParam,maxValue);

        double lower = minValue, upper = maxValue;
        double width = upper-lower;
        const double precision = width/10;

        while(width>precision)
        {
            cout << "\tbetween " << lower << " and " << upper << endl;
            Map[(2*lower+upper)/3]=getNbFromParam(image,config,nameParam,(2*lower+upper)/3);
            Map[(lower+2*upper)/3]=getNbFromParam(image,config,nameParam,(lower+2*upper)/3);
            if( abs(Map[(lower+2*upper)/3]-Map[lower]) < abs(Map[upper]-Map[(2*lower+upper)/3]) )
                upper = (2*lower+upper)/3;
            else
                lower = (lower+2*upper)/3;
            width = upper-lower;
        }
        //save the optimal value of the parameter
        //const float finalValue = (upper+lower)/2;
        //config.add(nameParam,finalValue);
    }
    else
    {
        for(double v=minValue;v<=maxValue;v+=stepValue)
        {
            Map[v]=getNbFromParam(image,config,nameParam,v);
        }
    }

    //display and save to file
    ofstream output_file((outputPath+"autoconfig_"+nameParam+".txt").c_str(), ios::out | ios::trunc);
    cout << "*\t" << nameParam << " results :\t*" <<endl;
    cout << nameParam << "\tN" << endl;
    output_file << nameParam << "\tN" << endl;
    for(map<double,int>::iterator i=Map.begin();i!=Map.end();++i)
    {
        cout << (*i).first << "\t" << (*i).second << endl;
        output_file << (*i).first << "\t" << (*i).second << endl;
    }
    output_file.close();

    //ofstream config_file("config.ini", ios::out | ios::trunc);
    //config_file << config;
    //config_file.close();
    cout << "done !" << endl;

    cout << endl << "**************************" << endl << endl;
}*/

    //}@
/*
int main()
{
    ConfigFile config("Tracker.ini");


    if(config.read("autoconfig",true))
    {
      autoconfig(config);
      return 0;
    }
    else
    {
        cout << "It is possible to configure automatically each parameter influencing the tracking. To do so, you have to add in the config.ini the following lines (exemple with the parameter minSep) : " << endl;
        cout << "\tautoconfig = 1" << endl;
        cout << "\tnameParam = minSep" << endl;
        cout << "\tminValue = 0.6" << endl;
        cout << "\tmaxValue = 0.85" << endl;
        cout << "The autoconfiguration routine will then vary the parameter between minValue and maxValue, looking for the region with least variation in the number of tracked particles.";
    }

    //read and parse the image name from the configuration file
    string Imagefile,file_prefix,file_sufix,outputPath, time_prefix;
    int nbdigit_z,nbdigit_t;
    if(!config.readInto<string>( Imagefile, "Imagefile" ))
      cerr << "the \"Imagefile\" entry is missing in config.ini";
    parseImageFileName(Imagefile,file_prefix,file_sufix,nbdigit_z,nbdigit_t);

    config.readInto<string>(outputPath,"outputPath","");

    //get the configuration data needed to open files
    int xsize, ysize,zsize,zoffset,tsize,toffset,channel;
    config.readInto(xsize,"xsize");
    config.readInto(ysize,"ysize");
    config.readInto(zsize,"zsize",1);
    config.readInto(zoffset,"zoffset",0);
    config.readInto(tsize,"tsize",1);
    config.readInto(toffset,"toffset",0);
    config.readInto(channel,"channel",0);

    float VoxelWidth,VoxelDepth,Zratio=1;
    if(zsize!=1 && (!config.readInto(VoxelWidth,"VoxelWidth") || !config.readInto(VoxelDepth,"VoxelDepth")) )
      cerr << "the \"VoxelWidth\" or \"VoxelDepth\" entry is missing in config.ini\nLook into the data from the microscope or XML export.";
    else
      Zratio = VoxelDepth / VoxelWidth;


    //Creates a 3D image of the size of each frame of the dataset
    CImg<float> image(xsize,ysize,zsize,1);


    const bool view = config.read("view",true);

    if(view)
    {
            cout << endl;
            cout << "* In a parameter calibration purpose, lots of images representing intermediate and final steps will be displayed." << endl;
            cout << "* To disable this calibration mode and for fast processing, edit the config.ini file in order to have the following line :" << endl;
            cout << "* view=0" << endl;
            cout << endl;
    }
    else
    {
            cout << endl;
            cout << "* If you want to display images representing intermediate and final steps, edit the config.ini file in order to have the following line :" << endl;
            cout << "* view=1" << endl;
            cout << endl;
    }


    if(nbdigit_t!=0)
      for(int t=0; t<tsize; ++t)
      {
        ostringstream oss;
        oss << t+toffset;
        string tstr =oss.str();

        for(int j = tstr.length();j<nbdigit_t;j++)
        {
         tstr = "0"+ tstr;
        }

        time_prefix = file_prefix + tstr;
        if(nbdigit_z!=0) time_prefix += "_z";
        cout << endl;
        cout << "\t\t\t*******************" << endl;
        cout << "\t\t\t***time " << tstr << " ***" << endl;
        cout << "\t\t\t*******************" << endl << endl;

        openRoutineZ(image,time_prefix,file_sufix,nbdigit_z,zoffset,channel);
        GraphicParticles Centers = trackXYZ(image,config);

        //export into a file
        Centers.exportToFile(outputPath+"t"+tstr+".dat");

      }
    else
    {
        openRoutineZ(image,file_prefix,file_sufix,nbdigit_z,zoffset,channel);
        //image.get_crop(0,image.dimy()/2,0,image.dimx()-1,image.dimy()/2,image.dimz()-1).permute_axes("xzyv").display().save("c:/bin/XZ.tif");

        GraphicParticles Centers = trackXYZ(image,config);

        //export into a file
        Centers.exportToFile(outputPath+"t0.dat");
    }


    return 0;
}*/
