#include "../tracker.hpp"

using namespace std;

#ifndef INSTAL_PATH
#define INSTAL_PATH "c:/bin/"
#endif

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout << "Syntax : aquireWisdom [path]filename" << endl;
		return EXIT_FAILURE;
	}
	string inputFile(argv[1]);
	try
	{
		auto_ptr<LifFile> lif(new LifFile(inputFile));
		const size_t serie = lif->chooseSerie();
		LifTracker track(lif,serie,lif->chooseChannel(serie),FFTW_EXHAUSTIVE);
		cout << "tracker ok"<<endl;
		LifTracker::saveWisdom(INSTAL_PATH "wisdom.fftw");
		cout<<"done !"<<endl;

	}
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
