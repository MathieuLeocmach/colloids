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


#ifndef lif_reader_H
#define lif_reader_H

#include "tinyxml.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>
#include <stdexcept>
#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace Colloids{
struct ChannelData;
struct DimensionData;
struct ScannerSettingRecord;
struct FilterSettingRecord;

class LifSerieHeader
{
    protected:
        std::string name;
        std::map<std::string, DimensionData> dimensions;
        std::vector<ChannelData> channels;
        std::vector<unsigned long long> timeStamps;
        std::map<std::string, ScannerSettingRecord> scannerSettings;
        TiXmlElement *rootElement;

    public:
        explicit LifSerieHeader(TiXmlElement *root);

        size_t chooseChannel() const;
        std::string getName() const {return this->name;};
        double getVoxelSize(const size_t d) const;
        int getResolution(const size_t channel) const;
        size_t getNbTimeSteps() const;
        std::vector<size_t> getSpatialDimensions() const;
        size_t getNbPixelsInOneTimeStep() const;
        size_t getNbPixelsInOneSlice() const;
        double getZXratio() const;
        const std::map<std::string, DimensionData>& getDimensionsData() const {return dimensions;};
        const std::vector<ChannelData>& getChannels() const {return channels;};

    private:
        void parseImage(TiXmlNode *elementImage);
        void parseImageDescription(TiXmlNode *elementImageDescription);
        void parseTimeStampList(TiXmlNode *elementTimeStampList);
        void parseHardwareSettingList(TiXmlNode *elementHardwareSettingList);
        void parseScannerSetting(TiXmlNode *elementScannerSetting);
};

class LifSerie : public LifSerieHeader, boost::noncopyable
{
    unsigned long long offset;
    unsigned long long memorySize;
    std::ifstream file;
    std::streampos fileSize;

    public:
        explicit LifSerie(LifSerieHeader serie, const std::string &filename, unsigned long long offset, unsigned long long memorySize);

        void fill3DBuffer(void* buffer, size_t t=0);
        void fill2DBuffer(void* buffer, size_t t=0, size_t z=0);
        std::istreambuf_iterator<char> begin(size_t t=0);
        std::streampos tellg(){return file.tellg();}

    private:
        unsigned long long getOffset(size_t t=0) const;
};

class LifHeader : boost::noncopyable
{
    int lifVersion;
    std::string name;

    protected:
        TiXmlDocument header;
        boost::ptr_vector<LifSerieHeader> series;

    public:
        explicit LifHeader(TiXmlDocument &header);
        explicit LifHeader(std::string &header);

        const TiXmlDocument& getXMLHeader() const{return this->header;};
        std::string getName() const {return this->name;};
        const int& getVersion() const {return this->lifVersion;};
        size_t getNbSeries() const {return this->series.size();}
        size_t chooseSerieNumber() const;
        LifSerieHeader& getSerieHeader(size_t s){return this->series[s];};
        const LifSerieHeader& getSerieHeader(size_t s) const {return this->series[s];};
        LifSerieHeader& chooseSerieHeader(){return getSerieHeader(chooseSerieNumber());};
        const LifSerieHeader& chooseSerieHeader() const {return getSerieHeader(chooseSerieNumber());};

    private:
        void parseHeader();
};

class LifReader : boost::noncopyable
{
    LifHeader *header;
    std::ifstream file;
    std::streampos fileSize;
    boost::ptr_vector<LifSerie> series;

    public:
        explicit LifReader(const std::string &filename);

        const LifHeader& getLifHeader() const {return *this->header;};
        const TiXmlDocument& getXMLHeader() const{return getLifHeader().getXMLHeader();};
        std::string getName() const {return getLifHeader().getName();};
        const int& getVersion() const {return getLifHeader().getVersion();};
        size_t getNbSeries() const {return this->series.size();}
        size_t chooseSerieNumber() const {return getLifHeader().chooseSerieNumber();};

        LifSerie& getSerie(size_t s){return series[s];};
        const LifSerie& getSerie(size_t s) const {return series[s];};
        LifSerie& chooseSerie(){return getSerie(chooseSerieNumber());};
        const LifSerie& chooseSerie() const {return getSerie(chooseSerieNumber());};

    private:
        int readInt();
        unsigned int readUnsignedInt();
        unsigned long long readUnsignedLongLong();


};

/** Struct of channel data*/
struct ChannelData
{
    int dataType; // 0 Integer, 1 Float
    int channelTag; // 0 Gray value, 1 Red, 2 Green, 3 Blue
    int resolution; // Bits per pixel
    std::string nameOfMeasuredQuantity;
    double minimum;
    double maximum;
    std::string unit;
    std::string LUTName;
    bool isLUTInverted; // 0 Normal LUT, 1 Inverted Order
    unsigned long long bytesInc; // Distance from the first channel in Bytes
    int bitInc;

    explicit ChannelData(TiXmlElement *element);
    inline const std::string getName() const
    {
        switch(channelTag)
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
    int dimID; // 0 Not valid, 1 X, 2 Y, 3 Z, 4 T, 5 Lambda, 6 Rotation, 7 XT Slices, 8 T Slices
    int numberOfElements; // Number of elements in this dimension
    double origin; // Physical position of the first element (left pixel side)
    double length; // Physical length from the first left pixel side to the last left pixel side
    std::string unit; // Physical unit
    unsigned long long bytesInc; // Distance from the one element to the next in this dimension
    int bitInc;

    explicit DimensionData(TiXmlElement *element);
    inline const std::string getName() const
    {
        switch(dimID)
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
    std::string identifier;
    std::string unit;
    std::string description;
    int data;
    std::string variant;
    int variantType;

    explicit ScannerSettingRecord(TiXmlElement *element);
};

/** Struct of Filter Setting */
struct FilterSettingRecord
{
    std::string objectName;
    std::string className;
    std::string attribute;
    std::string description;
    int data;
    std::string variant;
    int variantType;

    explicit FilterSettingRecord(TiXmlElement *element);
};

std::ostream& operator<< (std::ostream& out, const LifSerieHeader &s );
std::ostream& operator<< (std::ostream& out, const LifHeader &r );
std::ostream& operator<< (std::ostream& out, const DimensionData &r );
}
#endif
