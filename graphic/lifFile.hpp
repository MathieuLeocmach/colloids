/**
 * \file lifFile.cpp
 * \brief Define classes to read directly from Leica .lif file
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 31 March 2009
 *
 * Strongly inspired from BioimageXD VTK implementation but translated to fit with CImg framework
 *
 */


#ifndef lif_file_H
#define lif_file_H

#define TIXML_USE_STL 1
#include "../tinyxml/tinyxml.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>
#include <stdexcept>
#include <boost/utility.hpp>

/** Struct of channel data*/
struct ChannelData
{
    int DataType; // 0 Integer, 1 Float
    int ChannelTag; // 0 Gray value, 1 Red, 2 Green, 3 Blue
    int Resolution; // Bits per pixel
    const char *NameOfMeasuredQuantity;
    double Min;
    double Max;
    const char *Unit;
    const char *LUTName;
    int IsLUTInverted; // 0 Normal LUT, 1 Inverted Order
    unsigned long long BytesInc; // Distance from the first channel in Bytes
    int BitInc;

    explicit ChannelData(TiXmlElement *Element);
    inline const char * getName()
    {
        switch(ChannelTag)
        {
            case 0 : return "Gray";
            case 1 : return "Red";
            case 2 : return "Green";
            case 3 : return "Blue";
            default : return "Unknown color";
        }
    }
};

/** Struct of dimension data*/
struct DimensionData
{
    int DimID; // 0 Not valid, 1 X, 2 Y, 3 Z, 4 T, 5 Lambda, 6 Rotation, 7 XT Slices, 8 T Slices
    size_t NumberOfElements; // Number of elements in this dimension
    double Origin; // Physical position of the first element (left pixel side)
    double Length; // Physical length from the first left pixel side to the last left pixel side
    const char *Unit; // Physical unit
    unsigned long long BytesInc; // Distance from the one element to the next in this dimension
    int BitInc;

    explicit DimensionData(TiXmlElement *Element);
    inline const char * getName()
    {
        switch(DimID)
        {
            case 1 : return "X";
            case 2 : return "Y";
            case 3 : return "Z";
            case 4 : return "T";
            case 5 : return "Lambda";
            case 6 : return "Rotation";
            case 7 : return "XT Slices";
            case 8 : return "TSlices";
            default : return "Invalid";
        }
    };
};

/**
    \brief Struct of scanner setting
    Maybe because the original Leica program was written in MS Visual Basic,
    the type of the field "Variant" has to be defined from a test on the VariantType field :
    3:long
    5:double
    8:char*
    11:bool
    17:char (or int)
    See http://en.wikipedia.org/wiki/Variant_type
*/
struct ScannerSettingRecord
{
    const char *Identifier;
    const char *Unit;
    const char *Description;
    int Data;
    const char *Variant;
    int VariantType;

    explicit ScannerSettingRecord(TiXmlElement *Element);
};

/** Struct of Filter Setting */
struct FilterSettingRecord
{
    const char *ObjectName;
    const char *ClassName;
    const char *Attribute;
    const char *Description;
    int Data;
    const char *Variant;
    int VariantType;

    explicit FilterSettingRecord(TiXmlElement *Element);
};

class LifFile : boost::noncopyable
{
    std::ifstream file;
    std::streampos FileSize;

    int readInt();
    unsigned int readUnsignedInt();
    unsigned long long readUnsignedLongLong();

    int ParseInfoHeader(TiXmlElement *rootElement, int root= 1);
    void ParseImage(TiXmlNode *elementImage);
    void ParseImageDescription(TiXmlNode *elementImageDescription);
    void ParseTimeStampList(TiXmlNode *elementTimeStampList);
    void ParseHardwareSettingList(TiXmlNode *elementHardwareSettingList);
    void ParseScannerSetting(TiXmlNode *elementScannerSetting);

    unsigned long long GetTimePointOffset(const size_t & serie, const size_t & timepoint);

    public:
        TiXmlDocument Header;
        int LifVersion;
        std::deque< std::deque<DimensionData*> * > Dimensions;
        std::deque< std::deque<ChannelData*> * > Channels;
        std::deque< std::deque<unsigned long long> > TimeStamps;
        std::deque<const char*> Names;
        std::deque< std::map<std::string,ScannerSettingRecord*> > ScannerSettings;
        std::vector<unsigned long long> Offsets;
        std::vector<unsigned long long> ImageSizes;

        explicit LifFile(const std::string &filename);

        size_t chooseSerie();
        size_t chooseChannel(const size_t &serie);

        double getVoxelSize(const size_t &serie, const size_t &d);
        int GetResolution(const size_t &serie, const size_t &channel);
        void fill3DBuffer(void* buffer,const size_t &serie,const size_t &frame);
        unsigned int getNbPixelsInOneFrame(const size_t &serie);
        size_t getNbFrames(const size_t &serie);

};

#endif
