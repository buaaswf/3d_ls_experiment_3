#pragma once
#include <vector>
#include"ThreeDim_LevelSet.h"

//#include "statistics.h"
#define  input2 "D:\\data\\clean\\polypseginputdata\\origin\\"
#define output "D:\\swfdata20140420res\\polyp\\"
#define inputt "D:\\swfdata20140420res\\polyp\\thickness\\"
using namespace std;
class Seed
{
public:
	int x;
	int y;
	int z;
public:
	Seed(int x,int y,int z)
	{
		this->x=x;
		this->y=y;
		this->z=z;
	}
};
class Polyp
{
private:
	Raw *colon;
public:
		vector<Seed> seedlist;
		Polyp(void);
		~Polyp(void);
		vector<Seed> readSeedFromTXT(char *txtfilepath);
		Raw* initialRegion(vector<Seed> seedlist,int size,int l, int m, int n);
		void polypDetect(string dir);

		/**
		 \brief	Polyp detect.
		
		 \param	dir				The dir.
		 \param	dirthickness	The dirthickness.
		 \param	method			The method.
		 method==1 minimal surface
		 method==2 region growing
		 */

		void polypDetect(string dir,string dirthickness,string dirseg,int method,int i);
		void computePolyp(Raw *polyp,char* dir);

};

