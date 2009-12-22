/**
 * \file lifFile.cpp
 * \brief Implement classes to read directly from Leica .lif file
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 31 March 2009
 *
 * Strongly inspired from BioimageXD VTK implementation but translated to fit with CImg framework
 *
 */

#include <stdexcept>
#include "lifFile.hpp"

using namespace std;

/** \brief Constructor from lif file name */
LifFile::LifFile(const string &filename)
{
    //FILE * output;
    //output = fopen((filename+"_header.xml").c_str(),"w");
    //ofstream output((filename+"_header.xml").c_str(), ios::out | ios::trunc);
    file.open(filename.c_str(), ios::in | ios::binary);
    if(!file)
        throw invalid_argument(("No such file as "+filename).c_str());

    // Get size of the file
    file.seekg(0,ios::end);
    FileSize = file.tellg();
    file.seekg(0,ios::beg);

    const int MemBlockCode = 0x70,TestCode = 0x2a;
    if(readInt() != MemBlockCode)
        throw invalid_argument((filename+" is not a Leica SP5 file").c_str());

    // Skip the size of next block
    file.seekg(4,ios::cur);

    char lifChar;
    file>>lifChar;
    if(lifChar != TestCode)
        throw invalid_argument((filename+" is not a Leica SP5 file").c_str());

    unsigned int xmlChars = readUnsignedInt();
    //cout << xmlChars<<endl;
    xmlChars*=2;

    // Read and parse xml header
    string xmlString;
    xmlString.reserve(xmlChars/2);
	char *xmlHeader = new char[xmlChars];
    file.read(xmlHeader,xmlChars);
    for(unsigned int p=0;p<xmlChars/2;++p)
        xmlString[p] = xmlHeader[2*p];
    delete[] xmlHeader;
    Header.Parse(xmlString.c_str(),0);
    Header.RootElement()->QueryIntAttribute("Version",&LifVersion);
    //Header.Print(output);
    ParseInfoHeader(Header.RootElement());

    // Find image offsets
    Offsets.assign(Names.size(),0);
    ImageSizes.assign(Names.size(),0);
    size_t offsetId = 0;

    while (file.tellg() < FileSize)
    {
        // Check LIF test value
        int lifCheck = readInt();

        if (lifCheck != MemBlockCode) {
            long long fileOffset = file.tellg();
            fileOffset -= 4;
            cerr<< "ReadLIFHeader: File contains wrong MemBlockCode: " << lifCheck << " at: " << fileOffset<<endl;
            exit(1);
        }

        // Don't care about the size of the next block
        file.seekg(4,ios::cur);
        // Read testcode
        file>>lifChar;
        if (lifChar != TestCode)
        {
            long long fileOffset =file.tellg();
            fileOffset -= 1;
            cerr<< "ReadLIFHeader: File contains wrong TestCode: " << lifChar << " at: " << fileOffset<<endl;
            exit(1);
        }

        // Read size of memory, this is 4 bytes in version 1 and 8 in version 2
        unsigned long long memorySize;
        if (this->LifVersion >= 2)
            memorySize = readUnsignedLongLong();
        else
            memorySize = readUnsignedInt();


        // Find next testcode
        lifChar=0;
        while (lifChar != TestCode) {file>>lifChar;}
        unsigned int memDescrSize = readUnsignedInt() * 2;
        // Skip over memory description
        file.seekg(memDescrSize,ios::cur);

        // Add image offset if memory size is > 0
        if (memorySize > 0)
        {
            Offsets[offsetId]=file.tellg();
            ImageSizes[offsetId]=memorySize;
            offsetId++;
            file.seekg(static_cast<streampos>(memorySize),ios::cur);
        }
    }

    //output.close();
    //fclose(output);
    return;
}

/** \brief read an int form file advancing the cursor*/
int LifFile::readInt()
{
    char buffer[4];
    file.read(buffer,4);
    return *((int*)(buffer));
}
/** \brief read an unsigned int form file advancing the cursor*/
unsigned int LifFile::readUnsignedInt()
{
    char buffer[4];
    file.read(buffer,4);
    return *((unsigned int*)(buffer));
}
/** \brief read an unsigned long long int form file advancing the cursor*/
unsigned long long LifFile::readUnsignedLongLong()
{
    char buffer[8];
    file.read(buffer,8);
    return *((unsigned long long*)(buffer));
}

/** \brief Transforms XML Header into exploitable informations*/
int LifFile::ParseInfoHeader(TiXmlElement *rootElement, int root)
{
    TiXmlElement *elementElement = rootElement;
    if (root)
	{
        elementElement = rootElement->FirstChildElement();
        if (!elementElement)
        {
            cerr<< "ParseXMLHeader: No element Element found after root element."<<endl;
            return 0;
        }
	}
    TiXmlNode *elementData = elementElement->FirstChild("Data");//->ToElement();
    TiXmlNode *elementImage = elementData->FirstChild("Image");
    TiXmlNode *elementMemory = elementElement->FirstChild("Memory");
    TiXmlNode *elementChildren = elementElement->FirstChild("Children");


    // If Image element found
    if (elementImage)
    {
        this->ParseImage(elementImage);
        // Check that image info is read correctly and then add image name
        if (Channels.size() > Names.size() || Dimensions.size() > Names.size())
        {
            const char *imageName = elementElement->Attribute("Name");
            Names.push_back(imageName);
        }
    }

    // If Children element found
    if (elementChildren && !elementChildren->NoChildren())
    {
        TiXmlNode * child =0;
        while( child = elementChildren->IterateChildren("Element", child ) )
            this->ParseInfoHeader(child->ToElement(),0);
    }

    return 1;
}

/** \brief Transforms image XML header into exploitable informations*/
void LifFile::ParseImage(TiXmlNode *elementImage)
{
    TiXmlNode *elementImageDescription = elementImage->FirstChild("ImageDescription");
    TiXmlNode *elementTimeStampList = elementImage->FirstChild("TimeStampList");
    TiXmlNode *elementHardwareSettingList = 0;

    if (elementImageDescription)
        ParseImageDescription(elementImageDescription);

    // Parse time stamps even if there aren't any, then add empty
    // Unsigned Long Long to timestamps vector
    ParseTimeStampList(elementTimeStampList);

    TiXmlNode *Iterator = 0;
    string Name;
    // Parse Hardware Setting List even if there aren't
    while( Iterator = elementImage->IterateChildren("Attachment",Iterator ) )
    {
        Name = Iterator->ToElement()->Attribute("Name");
        if(Name=="HardwareSettingList" && !Iterator->NoChildren())
            elementHardwareSettingList = Iterator;
    }
    ParseHardwareSettingList(elementHardwareSettingList);
}

/** \brief Transforms XML describing an image into exploitable informations*/
void LifFile::ParseImageDescription(TiXmlNode *elementImageDescription)
{
    TiXmlNode *elementChannels = elementImageDescription->FirstChild("Channels");
    TiXmlNode *elementDimensions = elementImageDescription->FirstChild("Dimensions");
    if (elementChannels && elementDimensions)
	{
        deque<ChannelData*> *ImgChannels = new deque<ChannelData*>;
        deque<DimensionData*> *ImgDimensions = new deque<DimensionData*>;
        TiXmlNode *Iterator = 0;

        // Get info of channels
        while( Iterator = elementChannels->IterateChildren(Iterator ) )
            ImgChannels->push_back(new ChannelData(Iterator->ToElement()));

        // Get info of dimensions
        Iterator=0;
        while( Iterator = elementDimensions->IterateChildren(Iterator ) )
            ImgDimensions->push_back(new DimensionData(Iterator->ToElement()));

        Channels.push_back(ImgChannels);
        Dimensions.push_back(ImgDimensions);
	}
}
/** \brief Transforms XML describing time stamps into exploitable informations*/
void LifFile::ParseTimeStampList(TiXmlNode *elementTimeStampList)
{
    deque<unsigned long long> timeStampArray;
    if (elementTimeStampList)
    {
        unsigned long highInt;
        unsigned long lowInt;
        unsigned long long timeStamp;
        TiXmlNode *Iterator = 0;
        while( Iterator = elementTimeStampList->IterateChildren(Iterator ) )
        {
            Iterator->ToElement()->QueryValueAttribute<unsigned long>("HighInteger",&highInt);
            Iterator->ToElement()->QueryValueAttribute<unsigned long>("LowInteger",&lowInt);
            timeStamp = highInt;
            timeStamp <<= 32;
            timeStamp += lowInt;
            timeStamp /= 10000; // Convert to ms
            timeStampArray.push_back(timeStamp);
        }
    }
    TimeStamps.push_back(timeStampArray);
}
/** \brief Transforms XML describing Hardware settings into exploitable informations*/
void LifFile::ParseHardwareSettingList(TiXmlNode *elementHardwareSettingList)
{
    TiXmlNode *elementHardwareSetting = elementHardwareSettingList->FirstChild();
    TiXmlNode *elementScannerSetting = elementHardwareSetting->FirstChild("ScannerSetting");
    TiXmlNode *elementFilterSetting = elementHardwareSetting->FirstChild("FilterSetting");

    ParseScannerSetting(elementScannerSetting);
    //no real need of filter settings for the moment
    //ParseFilterSetting(elementFilterSetting);
}
/** \brief Transforms XML describing Scanner settings into exploitable informations*/
void LifFile::ParseScannerSetting(TiXmlNode *elementScannerSetting)
{
    map<string,ScannerSettingRecord*> ScannerSettingMap;
    if(elementScannerSetting)
    {
        ScannerSettingRecord *record;
        TiXmlNode *Iterator = 0;
        while( Iterator = elementScannerSetting->IterateChildren("ScannerSettingRecord",Iterator ) )
        {
            record = new ScannerSettingRecord(Iterator->ToElement());
            ScannerSettingMap.insert(pair<string,ScannerSettingRecord*>(record->Identifier,record));
        }
    }
    ScannerSettings.push_back(ScannerSettingMap);
}

/** \brief constructor from XML  */
ChannelData::ChannelData(TiXmlElement *Element)
{
    Element->QueryIntAttribute("DataType",&DataType);
    Element->QueryIntAttribute("ChannelTag",&ChannelTag);
    Element->QueryIntAttribute("Resolution",&Resolution);
    NameOfMeasuredQuantity = Element->Attribute("NameOfMeasuredQuantity");
    if (!NameOfMeasuredQuantity) NameOfMeasuredQuantity = "";
    Element->QueryDoubleAttribute("Min",&Min);
    Element->QueryDoubleAttribute("Max",&Max);
    Unit = Element->Attribute("Unit");
    if (!Unit) Unit = "";
    LUTName = Element->Attribute("LUTName");
    if (!LUTName) LUTName = "";
    Element->QueryIntAttribute("IsLUTInverted",&IsLUTInverted);
    // Bytes Inc is 64 bits in LIF version 2 but GetScalarAttribute only allows
    // maximum of unsigned long which can be 32 bits.
    Element->QueryValueAttribute<unsigned long long>("BytesInc",&BytesInc);
    Element->QueryIntAttribute("BitInc",&BitInc);
}

/** \brief constructor from XML  */
DimensionData::DimensionData(TiXmlElement *Element)
{
    Element->QueryIntAttribute("DimID",&DimID);
    Element->QueryIntAttribute("NumberOfElements",(int*)&NumberOfElements);
    Element->QueryDoubleAttribute("Origin",&Origin);
    Element->QueryDoubleAttribute("Length",&Length);
    Unit = Element->Attribute("Unit");
    if (!Unit) Unit = "";
    // Bytes Inc is 64 bits in LIF version 2 but GetScalarAttribute only allows
    // maximum of unsigned long which can be 32 bits.
    Element->QueryValueAttribute<unsigned long long>("BytesInc",&BytesInc);
    Element->QueryIntAttribute("BitInc",&BitInc);
}

/** \brief constructor from XML  */
ScannerSettingRecord::ScannerSettingRecord(TiXmlElement *Element)
{
    Identifier = Element->Attribute("Identifier");
    Unit = Element->Attribute("Unit");
    if (!Unit) Unit = "";
    Description = Element->Attribute("Description");
    Element->QueryIntAttribute("Data",&Data);
    Variant = Element->Attribute("Variant");
    Element->QueryIntAttribute("VariantType",&VariantType);
}

/** \brief constructor from XML  */
FilterSettingRecord::FilterSettingRecord(TiXmlElement *Element)
{
    ObjectName = Element->Attribute("ObjectName");
    ClassName = Element->Attribute("ClassName");
    Attribute = Element->Attribute("Attribute");
    Description = Element->Attribute("Description");
    Element->QueryIntAttribute("Data",&Data);
    Variant = Element->Attribute("Variant");
    Element->QueryIntAttribute("VariantType",&VariantType);
}
/** \brief get the real size of a pixel (in meters) in the dimension d, for the serie serie */
double LifFile::getVoxelSize(const size_t &serie, const size_t & d)
{
    const string Voxel = "dblVoxel", dims[3] = {"X","Y","Z"};
    stringstream sstream( ScannerSettings[serie][Voxel+dims[d]]->Variant );
    double VoxelSize;
    sstream >> VoxelSize;
    return VoxelSize;
}

/** \brief get the resolution of a channel of a serie, in Bits per pixel*/
int LifFile::GetResolution(const size_t & serie, const size_t & channel)
{
    return Channels[serie]->at(channel)->Resolution;
}

/** \brief get the offset between the begining of the serie and the begining of the frame*/
unsigned long long LifFile::GetTimePointOffset(const size_t & serie, const size_t & timepoint)
{
  for (deque<DimensionData*>::const_iterator dimIter = Dimensions[serie]->begin();
       dimIter != Dimensions[serie]->end(); dimIter++)
    {
      if ((*dimIter)->DimID == 4) return (*dimIter)->BytesInc * timepoint;
    }

  return 0;
}

unsigned int LifFile::getNbPixelsInOneFrame(const size_t & serie)
{
    if (serie >= Dimensions.size()) return 0;

    unsigned int voxels = 1;
    for (deque<DimensionData*>::const_iterator dimIter = Dimensions[serie]->begin();
        dimIter != Dimensions[serie]->end(); dimIter++)
    {
        if ((*dimIter)->DimID <4) voxels *= (*dimIter)->NumberOfElements;
    }
    return voxels;
}

size_t LifFile::getNbFrames(const size_t &serie)
{
    if (serie >= Dimensions.size()) return 0;
    if (Dimensions[serie]->size()<4) return 1;
    return Dimensions[serie]->at(3)->NumberOfElements;
}

/**
    \brief fill a memory buffer that already has the good dimension
    If more than one channel, all channels are retrieved interlaced

    */
void LifFile::fill3DBuffer(void* buffer,const size_t & serie, const size_t & frame)
{
    char *pos = static_cast<char*>(buffer);
    unsigned long long serieOffset = Offsets[serie];
    serieOffset += GetTimePointOffset(serie,frame);
    unsigned long int frameDataSize = getNbPixelsInOneFrame(serie)*Channels[serie]->size();
    file.seekg(serieOffset,ios::beg);
    file.read(pos,frameDataSize);
}

/** \brief Output the list of the series in the LIF file and ask the user to choose one */
size_t LifFile::chooseSerie()
{
    for(size_t i=0;i<Names.size();++i)
    {
        cout<<"("<<i<<")" << Names[i]<<":\t"
            <<Dimensions[i]->size()<<" dimensions\t";
        for(size_t d=0;d<Dimensions[i]->size();++d)
            cout<<Dimensions[i]->at(d)->getName()<<" "
                <<Dimensions[i]->at(d)->NumberOfElements<<", ";
        cout<<endl;
        cout<<"\tVoxel sizes (nm)\t";
        for(size_t d=0;d<3 && d<Dimensions[i]->size() ;++d)
            if(Dimensions[i]->at(d)->DimID < 4)
                cout<<getVoxelSize(i,d)*1E9<<"\t";
        cout<<endl;
        for(size_t ch=0;ch<Channels[i]->size();++ch)
            cout << "\t"<<Channels[i]->at(ch)->getName()<<endl;
    }
    size_t serie;
    do
    {
        cout<<"Chose a serie between 0 and "<<Names.size()-1<<": ";
        cin>>serie;
    }while(serie>=Names.size());
    return serie;
}

/** \brief allows the user to choose a channel */
size_t LifFile::chooseChannel(const size_t &serie)
{
    if(Channels[serie]->size()==1) return 0;
    size_t channel;
    for(size_t ch=0;ch<Channels[serie]->size();++ch)
        cout << "("<<ch<<")"<<Channels[serie]->at(ch)->getName()<<endl;
    do
    {
        cout<<"Chose a channel between 0 and "<<Channels[serie]->size()-1<<": ";
        cin>>channel;
    }while(channel>=Channels[serie]->size());
    return channel;
}

