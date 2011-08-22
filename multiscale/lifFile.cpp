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

 * \file lifFile.cpp
 * \brief Define classes to read directly from Leica .lif file
 * \author Mathieu Leocmach
 * \date 28 March 2009
 *
 * Strongly inspired from BioimageXD VTK implementation but using standard library instead of VTK objects
 *
 */
#include "lifFile.hpp"
#include <algorithm>
#include <numeric>
#include <sstream>

using namespace std;
using namespace Colloids;

/** @brief LifSerieHeader constructor  */
LifSerieHeader::LifSerieHeader(TiXmlElement *root) : name(root->Attribute("Name")), rootElement(root)
{
    //XML parsing
    TiXmlNode *elementImage = rootElement->FirstChild("Data")->FirstChild("Image");
    // If Image element found
    if (elementImage)
        parseImage(elementImage);
    return;
}

/** @brief Let the user choose a Channel  */
size_t LifSerieHeader::chooseChannel() const
{
    for(size_t ch=0;ch<channels.size();++ch)
        cout << "("<<ch<<")"<<channels[ch].getName()<<endl;

	if(channels.size() == 1)
		return 0;

	size_t channel;
    do
    {
        cout<<"Chose a channel between 0 and "<<channels.size()-1<<": ";
        cin>>channel;
    }while(channel >= channels.size());
    return channel;
}

/** \brief get the real size of a pixel (in meters) in the dimension d */
double LifSerieHeader::getVoxelSize(const size_t d) const
{
    const string voxel = "dblVoxel", dims = "XYZ";
    map<string, ScannerSettingRecord>::const_iterator it = scannerSettings.find(voxel+dims[d]);
    if(it == scannerSettings.end())
        return 0.0;
    stringstream sstream( it->second.variant);
    double voxelSize;
    sstream >> voxelSize;
    return voxelSize;
}

/** \brief get the resolution of a channel of a serie, in Bits per pixel*/
int LifSerieHeader::getResolution(const size_t channel) const
{
    return channels[channel].resolution;
}

/** @brief get the number of time steps  */
size_t LifSerieHeader::getNbTimeSteps() const
{
    map<string, DimensionData>::const_iterator it = dimensions.find("T");
    if(it == dimensions.end()) return 1;
    return it->second.numberOfElements;
}

/** @brief getSpatialDimensions  */
vector<size_t> LifSerieHeader::getSpatialDimensions() const
{
    vector<size_t> dims;
    for(map<string, DimensionData>::const_iterator d = dimensions.begin(); d != dimensions.end(); d++)
        if(d->second.dimID <4)
			dims.push_back(d->second.numberOfElements);
    return dims;
}

/** @brief get the number of pixels in one time step, ie the product of the spatial dimensions  */
size_t LifSerieHeader::getNbPixelsInOneTimeStep() const
{
    const vector<size_t> dims = getSpatialDimensions();
    return accumulate(dims.begin(), dims.end(), 1, multiplies<size_t>());
}

/** @brief get the number of pixels in one slice, ie the product of the two first spatial dimensions  */
size_t LifSerieHeader::getNbPixelsInOneSlice() const
{
    vector<size_t> dims = getSpatialDimensions();
    dims.resize(2, 1);
    return accumulate(dims.begin(), dims.end(), 1, multiplies<size_t>());
}

/** @brief get the ratio (Z voxel size)/(X voxel size)  */
double LifSerieHeader::getZXratio() const
{
    double x = getVoxelSize(0), z = getVoxelSize(2);
    if(z==0.0)
        return 1.0;
    else
        return z/x;
}



/** @brief parse the "Image" node of the XML header  */
void LifSerieHeader::parseImage(TiXmlNode *elementImage)
{
    TiXmlNode *elementImageDescription = elementImage->FirstChild("ImageDescription");
    if (elementImageDescription)
        parseImageDescription(elementImageDescription);


    // Parse time stamps even if there aren't any, then add empty
    // Unsigned Long Long to timestamps vector
    TiXmlNode *elementTimeStampList = elementImage->FirstChild("TimeStampList");
    parseTimeStampList(elementTimeStampList);


    // Parse Hardware Setting List even if there aren't
    TiXmlNode *elementHardwareSettingList = 0;
    TiXmlNode *child = 0;
    string Name;
    while( (child = elementImage->IterateChildren("Attachment", child)) )
    {
        Name = child->ToElement()->Attribute("Name");
        if(Name=="HardwareSettingList" && !child->NoChildren())
            elementHardwareSettingList = child;
    }
    parseHardwareSettingList(elementHardwareSettingList);
}

/** @brief parse the "ImageDescription" node of the XML header  */
void LifSerieHeader::parseImageDescription(TiXmlNode *elementImageDescription)
{
    TiXmlNode *elementChannels = elementImageDescription->FirstChild("Channels");
    TiXmlNode *elementDimensions = elementImageDescription->FirstChild("Dimensions");
    if (elementChannels && elementDimensions)
	{
        TiXmlNode *child = 0;

        // Get info of channels
        while((child = elementChannels->IterateChildren(child)))
            this->channels.push_back(ChannelData(child->ToElement()));

        // Get info of dimensions
        child = 0;
        while((child = elementDimensions->IterateChildren(child)))
        {
            DimensionData data(child->ToElement());
            this->dimensions.insert(
                make_pair(data.getName(),data)
                );
        }
	}
}

/** @brief parse the "TimeStampList" node of the XML header  */
void LifSerieHeader::parseTimeStampList(TiXmlNode *elementTimeStampList)
{
    if (elementTimeStampList)
    {
        unsigned long highInt;
        unsigned long lowInt;
        unsigned long long timeStamp;
        TiXmlNode *child = 0;
        while((child = elementTimeStampList->IterateChildren(child)))
        {
            child->ToElement()->QueryValueAttribute<unsigned long>("HighInteger",&highInt);
            child->ToElement()->QueryValueAttribute<unsigned long>("LowInteger",&lowInt);
            timeStamp = highInt;
            timeStamp <<= 32;
            timeStamp += lowInt;
            timeStamp /= 10000; // Convert to ms
            this->timeStamps.push_back(timeStamp);
        }
    }
}

/** @brief parse the "HardwareSettingList" node of the XML header  */
void LifSerieHeader::parseHardwareSettingList(TiXmlNode *elementHardwareSettingList)
{
    TiXmlNode *elementHardwareSetting = elementHardwareSettingList->FirstChild();
    TiXmlNode *elementScannerSetting = elementHardwareSetting->FirstChild("ScannerSetting");

    parseScannerSetting(elementScannerSetting);
    //no real need of filter settings for the moment
    //TiXmlNode *elementFilterSetting = elementHardwareSetting->FirstChild("FilterSetting");
    //ParseFilterSetting(elementFilterSetting);
}

/** @brief parse the "ScannerSetting" node of the XML header  */
void LifSerieHeader::parseScannerSetting(TiXmlNode *elementScannerSetting)
{
    if(elementScannerSetting)
    {
        TiXmlNode *child = 0;
        while((child = elementScannerSetting->IterateChildren("ScannerSettingRecord", child)))
            this->scannerSettings.insert(
                make_pair(
                    child->ToElement()->Attribute("Identifier"),
                    ScannerSettingRecord(child->ToElement())
                    )
                );
    }
}

/** @brief LifSerie constructor  */
LifSerie::LifSerie(LifSerieHeader serie, const std::string &filename, unsigned long long offset, unsigned long long memorySize) :
 LifSerieHeader(serie)
{
    file.open(filename.c_str(), ios::in | ios::binary);
    if(!file)
        throw invalid_argument(("No such file as "+filename).c_str());

    //get size of the file
    file.seekg(0,ios::end);
    fileSize = file.tellg();
    //cout<<"offset: "<<offset<<"\tmemorySize: "<<memorySize<<"\tfileSize: "<<fileSize<<endl;
    //check the validity of the offset and memorysize parameters
    if(offset >= (unsigned long long)fileSize)
        throw invalid_argument("Offset is larger than file size");
    if(offset+memorySize > (unsigned long long)fileSize)
        throw invalid_argument("The end of the serie is further than the end of file");

    this->offset = offset;
    this->memorySize = memorySize;
}

/**
    \brief fill a memory buffer that already has the good dimension
    If more than one channel, all channels are retrieved interlaced
  */
void LifSerie::fill3DBuffer(void* buffer, size_t t)
{
    char *pos = static_cast<char*>(buffer);
    unsigned long int frameDataSize = getNbPixelsInOneTimeStep()*channels.size();
    file.seekg(getOffset(t) ,ios::beg);
    file.read(pos,frameDataSize);
}

/**
    \brief fill a memory buffer that already has the good dimension
    If more than one channel, all channels are retrieved interlaced
  */
void LifSerie::fill2DBuffer(void* buffer, size_t t, size_t z)
{
    char *pos = static_cast<char*>(buffer);
    unsigned long int sliceDataSize = getNbPixelsInOneSlice()*channels.size();
    file.seekg(getOffset(t) + z *  sliceDataSize, ios::beg);
    file.read(pos, sliceDataSize);
}

/** @brief return an iterator to the begining of the data of time step t
    No gestion of multi-channel.
*/
istreambuf_iterator<char> LifSerie::begin(size_t t)
{
    file.seekg(getOffset(t));
    return istreambuf_iterator<char>(file);
}


/** @brief get the position in file where starts the data of time step t  */
unsigned long long LifSerie::getOffset(size_t t) const
{
    if(t >= getNbTimeSteps())
        throw out_of_range("Time step out of range");
    if(t==0)
        return this->offset;

    map<string, DimensionData>::const_iterator it = dimensions.find("T");
    if(it == dimensions.end())
        return this->offset;

    return this->offset + it->second.bytesInc * t;
}

/** @brief LifHeader constructor  */
 LifHeader::LifHeader(TiXmlDocument &header)
{
    this->header = header;
    parseHeader();
    return;
}

/** @brief LifHeader
  *
  * @todo: document this function
  */
 LifHeader::LifHeader(std::string &header)
{
    this->header.Parse(header.c_str(),0);
    parseHeader();
    return;
}



/** @brief Display the available series and let the user choose one  */
size_t LifHeader::chooseSerieNumber() const
{
    cout<< *this <<endl;
    if(series.size() == 1)
		return 0;

	size_t s = 0;
    do
    {
        cout<<"Chose a serie between 0 and "<<series.size()-1<<": ";
        cin>>s;
    }while(s >= series.size());
    return s;
}

/** @brief parse the XML header  */
void LifHeader::parseHeader()
{
    header.RootElement()->QueryIntAttribute("Version", &lifVersion);

    TiXmlElement *elementElement = header.RootElement()->FirstChildElement();
    if (!elementElement)
        throw logic_error("No element Element found after root element.");

    this->name = elementElement->Attribute("Name");
    TiXmlNode *elementChildren = elementElement->FirstChild("Children");
    if (!elementChildren || elementChildren->NoChildren())
        throw logic_error("No serie in this experiment.");

    TiXmlNode *serieNode = 0;
    while( (serieNode = elementChildren->IterateChildren("Element", serieNode)))
        series.push_back(new LifSerieHeader(serieNode->ToElement()));
}



/** \brief Constructor from lif file name */
LifReader::LifReader(const string &filename) : file(filename.c_str(), ios::in | ios::binary)
{
    const int MemBlockCode = 0x70, TestCode = 0x2a;
    char lifChar;

    if(!file.is_open())
        throw invalid_argument(("No such file as "+filename).c_str());
    // Get size of the file
    file.seekg(0,ios::end);
    fileSize = file.tellg();
    file.seekg(0,ios::beg);
    cout<<filename<<" is "<<fileSize<<" bytes"<<endl;

    if(readInt() != MemBlockCode)
        throw invalid_argument((filename+" is not a Leica SP5 file").c_str());

    // Skip the size of next block
    file.seekg(4,ios::cur);
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
        xmlString.push_back(xmlHeader[2*p]);
    delete[] xmlHeader;

    header = new LifHeader(xmlString);

    size_t s = 0;
    while (file.tellg() < fileSize)
    {
        // Check LIF test value
        int lifCheck = readInt();
        if (lifCheck != MemBlockCode)
            throw logic_error("File contains wrong MemBlockCode");

        // Don't care about the size of the next block
        file.seekg(4,ios::cur);
        // Read testcode
        file>>lifChar;
        if (lifChar != TestCode)
            throw logic_error("File contains wrong TestCode");

        // Read size of memory, this is 4 bytes in version 1 and 8 in version 2
        unsigned long long memorySize;
        if (this->getVersion() >= 2)
            memorySize = readUnsignedLongLong();
        else
            memorySize = readUnsignedInt();


        // Find next testcode
        lifChar=0;
        while (lifChar != TestCode) {file>>lifChar;}
        unsigned int memDescrSize = readUnsignedInt() * 2;
        // Skip over memory description
        file.seekg(memDescrSize,ios::cur);

        // Add serie if memory size is > 0
        if (memorySize > 0)
        {
            if(s >= getLifHeader().getNbSeries())
                throw logic_error("Too many memory blocks");

            series.push_back(new LifSerie(
                    this->header->getSerieHeader(s),
                    filename,
                    file.tellg(),
                    memorySize)
                    );
            s++;
            //jump to the next memory block
            file.seekg(static_cast<streampos>(memorySize),ios::cur);
        }
    }

    return;
}

/** \brief read an int form file advancing the cursor*/
int LifReader::readInt()
{
    char buffer[4];
    file.read(buffer,4);
    return *((int*)(buffer));
}
/** \brief read an unsigned int form file advancing the cursor*/
unsigned int LifReader::readUnsignedInt()
{
    char buffer[4];
    file.read(buffer,4);
    return *((unsigned int*)(buffer));
}
/** \brief read an unsigned long long int form file advancing the cursor*/
unsigned long long LifReader::readUnsignedLongLong()
{
    char buffer[8];
    file.read(buffer,8);
    return *((unsigned long long*)(buffer));
}


/** \brief constructor from XML  */
ChannelData::ChannelData(TiXmlElement *element)
{
    element->QueryIntAttribute("DataType",&dataType);
    element->QueryIntAttribute("ChannelTag",&channelTag);
    element->QueryIntAttribute("Resolution",&resolution);
    nameOfMeasuredQuantity = element->Attribute("NameOfMeasuredQuantity");
    element->QueryDoubleAttribute("Min",&minimum);
    element->QueryDoubleAttribute("Max",&maximum);
    unit = element->Attribute("Unit");
    LUTName = element->Attribute("LUTName");
    element->QueryIntAttribute("IsLUTInverted",(int*)&isLUTInverted);
    // Bytes Inc is 64 bits in LIF version 2 but GetScalarAttribute only allows
    // maximum of unsigned long which can be 32 bits.
    element->QueryValueAttribute<unsigned long long>("BytesInc",&bytesInc);
    element->QueryIntAttribute("BitInc",&bitInc);
}

/** \brief constructor from XML  */
DimensionData::DimensionData(TiXmlElement *element)
{
    element->QueryIntAttribute("DimID",&dimID);
    element->QueryIntAttribute("NumberOfElements", &numberOfElements);
    element->QueryDoubleAttribute("Origin",&origin);
    element->QueryDoubleAttribute("Length",&length);
    unit = element->Attribute("Unit");
    // Bytes Inc is 64 bits in LIF version 2 but GetScalarAttribute only allows
    // maximum of unsigned long which can be 32 bits.
    element->QueryValueAttribute<unsigned long long>("BytesInc",&bytesInc);
    element->QueryIntAttribute("BitInc",&bitInc);
}

/** \brief constructor from XML  */
ScannerSettingRecord::ScannerSettingRecord(TiXmlElement *element)
{
    identifier = element->Attribute("Identifier");
    unit = element->Attribute("Unit");
    description = element->Attribute("Description");
    element->QueryIntAttribute("Data",&data);
    variant = element->Attribute("Variant");
    element->QueryIntAttribute("VariantType",&variantType);
}

/** \brief constructor from XML  */
FilterSettingRecord::FilterSettingRecord(TiXmlElement *element)
{
    objectName = element->Attribute("ObjectName");
    className = element->Attribute("ClassName");
    attribute = element->Attribute("Attribute");
    description = element->Attribute("Description");
    element->QueryIntAttribute("Data",&data);
    variant = element->Attribute("Variant");
    element->QueryIntAttribute("VariantType",&variantType);
}

/** @brief operator<<
  *
  * @todo: document this function
  */
ostream & Colloids::operator<<(ostream& out, const LifSerieHeader &s)
{
    out<< s.getName() <<":\t";
    const map<string, DimensionData> & dims = s.getDimensionsData();
    out<< dims.size()<<" dimensions\t";
    for(map<string, DimensionData>::const_iterator it=dims.begin(); it!=dims.end();++it)
        out<<it->second<<", ";
    out<<endl;
    out<<"\tVoxel sizes (nm)\t";
    for(map<string, DimensionData>::const_iterator it=dims.begin();it!=dims.end();++it)
        if(it->second.dimID < 4)
            cout<<s.getVoxelSize(it->second.dimID-1)*1E9<<"\t";
    out<<endl;
    const vector<ChannelData> &channels = s.getChannels();
    transform(
        channels.begin(), channels.end(),
        ostream_iterator<string>(out, "\t"),
        mem_fun_ref(&ChannelData::getName)
        );
    return out;
}

/** @brief operator<<
  *
  * @todo: document this function
  */
ostream & Colloids::operator<<(ostream& out, const LifHeader &r)
{
    out<< r.getName() <<endl;
    for(size_t s=0; s<r.getNbSeries(); ++s)
        out<<"("<<s<<") "<<r.getSerieHeader(s)<<endl;
    return out;
}

/** @brief operator<<
  *
  * @todo: document this function
  */
ostream & Colloids::operator<<(ostream& out, const DimensionData &d)
{
    out<<d.getName()<<" "<<d.numberOfElements;
    return out;
}


