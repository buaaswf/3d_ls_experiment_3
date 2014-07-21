#include <iostream>
#include <fstream>
#include <string>
#include "Polyp.h"
#include <queue>
//#import <c:\program files\common files\system\ado\msado15.dll>  
//
//rename("EOF", "adoEOF");//把"EOF"定义改为"adoEOF",避免和文件未eof冲突 	
//上面头文件的请照写,要不就出现很多错误 
//msado15.dll加载后会在我们c++工程项目里自动产生几个头文件,定义了下面 莫名其妙 的宏 类型 等等 
using namespace std;
void seedGrowing(Raw * src,vector<Seed> seedlist);
void seedGrowing(Raw * src, Raw *origin, vector<Seed> seedlist);
void seedGrowing(Raw * src, vector<Seed> seedlist, int threshold);
vector<Seed> seedlistdata(const int size);
Polyp::Polyp(void)
{
}


Polyp::~Polyp(void)
{
}
int CountLines(char *filename)
{
	ifstream ReadFile;
	int n=0;
	string tmp;
	ReadFile.open(filename,ios::in);//ios::in :readonly
	if(ReadFile.fail())// open fails
	{
		return 0;
	}
	else//if exit
	{
		while(getline(ReadFile,tmp,'\n'))
		{
			n++;
		}
		ReadFile.close();
		return n;
	}
}
vector<Seed> readSeedFromTXT(char *txtfilepath)
{
		ifstream file;
		int LINES;
		char filename[512]="inFile.txt";
		file.open(filename,ios::in);
		if(file.fail())
		{
			cout<<"not exit."<<endl;
			file.close();
		}
		else//if exit
		{
			LINES=CountLines(filename);
			int *tempInt=new int[LINES];
			char *tempChar=new char[LINES];
			int i=0;
			while(!file.eof()) //read data to array
			{

				file>>tempInt[i];
				file>>tempChar[i];
				i++;
			}
			file.close(); //close file
			for(i=0;i<LINES;i++)
				cout<<tempInt[i]<<"\t"<<tempChar[i]<<endl;
			delete []tempInt;
			delete []tempChar;
		}

		vector<Seed> seed;
		return seed;

}
Raw* Polyp::initialRegion(vector<Seed> seedlist,int size,int l, int m, int n)
{
	Raw *initialRegion=new Raw(l,m,n);
	//memset(initialRegion->getdata(),2,initialRegion->size());
	for (int i = 0; i < l*m*n; i++)
	{
		initialRegion->putXYZ(i,-2);
	}
	for (vector<Seed>::iterator it = seedlist.begin();it!=seedlist.end();++it)
	{
		
		//Seed *seed=new Seed(it->x,it->y,it->z);
		for (int i = it->x - size; i < it->x + size; ++i )
		{
			for (int j = it->y - size; j < it->y + size; ++j)
			{
				for (int k = it->z - size; k< it->z + size; ++k)
				{
					initialRegion->put(i,j,k,2);
				}
				
			}
			
		}
		
	}
	return initialRegion;
}
void Polyp::polypDetect(string dir)
{
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=20;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	unsigned char * thicknessdata=test.readStreamuchar(dirhead,l,m,n,212);
	n=50;
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	for (int i = 0; i < l*m*n; i++)
	{
		inputo[i]=(float) indata[i];		
	}

	Raw *input=new Raw(l,m,n,inputo);
	/*Filter *f=new Filter();*/
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//20140405 delete because of the existance of 
	ls->initialg(*input);
	//for (int i=0; i<input->getXsize(); i++)
	//{
	//	for (int j=0; j<input->getYsize(); j++)
	//	{
	//		for (int k=0; k<input->getZsize(); k++)
	//		{
	//			if (input->get(i,j,k)>=1)
	//			{
	//				initial->put(i,j,k,-2);
	//			}
	//			else initial->put(i,j,k,2);

	//		}
	//	}

	//}
	vector<Seed> seedlist;
	Seed *seed=new Seed(164, 373, 20);
	seedlist.push_back(*seed);
	Raw *initialdata=initialRegion(seedlist,20,l,m,n);
	*initial=ls->minimal_surface(*initialdata,*input,5.0,0.1,3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//Raw temp(*initial);
	//ls->outerwallauto(*initial,*input,10,0.1,-6,1.5,1,10,pt);
	//*initial -=temp;
	//char *outname2="outer5-8_2_20140405.raw";
	//char outdir2[200]=output;
	//strcat(outdir2,dirbody);
	//strcat(outdir2,outname2);
	//test.writeImageName(*initial,outdir2);
	//evaluate(dir,l,m,n);
}
vector<Seed> probabiltyouterpointleastcost(Raw *data, Seed *center,int range)
{
	int x = center->x;
	int y = center->y;
	int z = center->z;
	vector<Seed> datalist;
	for (int i = 0; i < data->getXsize();i++)
	{
		for (int j = 0; j < data->getYsize(); j++)
		{
			for (int k = 0; k < data->getZsize(); k++)
			{
				if ((i - x)*(i - x) + (j - y)*(j - y) + (k - z)*(k - z) <= range*range)
				{
					Seed *seed = new Seed(i, j, k);
					datalist.push_back(*seed);
				}
			}
		}
	}
	return datalist;


}
Seed * findLeastCostSeed(vector<Seed> list,Seed *center)
{
	int distance = 100000;
	int pos = 0;
	for (int i = 0; i < list.size(); i++)
	{
		int val = (list[i].x - center->x)*(list[i].x - center->x) + (list[i].y - center->y)*(list[i].y - center->y)
			+ (list[i].z - center->z)*(list[i].z - center->z);
		val < distance*distance ? distance = sqrt(val) : distance = distance;
		pos = i;
	}
	return &(list[pos]);
}
void Polyp::polypDetect(string dir,string dirthickness,string dirseg,int methodo,int i)
{
	
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=40;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\
	
	char dirheadt[200]=inputt; 
	char dirheadseg[200] = "D:\\data\\segmention\\polypseg\\";
	char dirbodyt[100];
	char dirbodyseg[100];
	char dirbody[100];
	//char dirbodyseg[100];
	strcpy(dirbody,dir.c_str());
	strcpy(dirbodyt,dirthickness.c_str());
	strcpy(dirbodyseg, dirseg.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	cout <<"dirbody thickness"<<dirbodyt<<endl;
	cout << "dirbody seg" <<dirbodyseg<<endl;
	strcat(dirhead,dirbody);
	strcat(dirheadt,dirbodyt);
	strcat(dirheadseg, dirbodyseg);
	cout << "dirhead" << dirhead<< endl;
	vector<Seed> list;
	list = seedlistdata(42);
	vector<Seed> seedlist;
	//Seed *seed = new Seed(248, 208, 15);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18:
	Seed *seed = new Seed(list[i]); //to be c
	float rate = 0;
	test.readStreamseg(dirheadseg, &rate);
	seed->z = seed->z*rate;
	int offset = seed->z - 20;
	short * indata = test.readStreamfseek(dirhead, &l, &m, &n, offset);

	seed->z = 20;
	
	unsigned char * thicknessdata=test.readStreamuchar(dirheadt,l,m,n,offset);//3036P 212 3033P 90 3033P peduncated:90 sessile1:70 sessile2:248,204,58:60
	//PIXTYPE * 
	n=50;
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];

	//for (int i=0; i<input->getXsize(); i++)
	//{
	//	for (int j=0; j<input->getYsize(); j++)
	//	{
	//		for (int k=0; k<input->getZsize(); k++)
	//		{
	//			if (input->get(i,j,k)>=1)
	//			{
	//				initial->put(i,j,k,-2);
	//			}
	//			else initial->put(i,j,k,2);

	//		}
	//	}

	//}

	seedlist.push_back(*seed);
	int method=methodo;
	if (method==1)
	{
		for (int i = 0; i < l*m*n; i++)
		{
			inputo[i] = (float)indata[i];
		}
		delete[] thicknessdata;
		Raw *input = new Raw(l, m, n, inputo);
		/*Filter *f=new Filter();*/
		//input=f->guass3DFilter(input,3);
		RawImage *write = new RawImage();
		ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
		//20140405 delete because of the existance of 
		ls->initialg(*input);
	Raw *initialdata=initialRegion(seedlist,2,l,m,n);
	*initial=ls->minimal_surface(*initialdata,*input,2*5.0,0.1,-3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//delete[] input;
	//Raw temp(*initial);
	//ls->outerwallauto(*initial,*input,10,0.1,-6,1.5,1,10,pt);
	//*initial -=temp;
	//char *outname2="outer5-8_2_20140405.raw";
	//char outdir2[200]=output;
	//strcat(outdir2,dirbody);
	//strcat(outdir2,outname2);
	//test.writeImageName(*initial,outdir2);
	//evaluate(dir,l,m,n);
	}
	else if (method==2)
	{
		for (int i = 0; i < l*m*n; i++)
		{
			inputo[i] = (float)thicknessdata[i];
		}
		PIXTYPE *origin = new PIXTYPE[l*m*(n+offset)];
		for (int j = 0; j < l*m*(n+offset); j++)
		{
			origin[j] = (float)indata[j];
		}

		Raw *origindata = new Raw(l, m, offset+n, origin);
		//delete[] thicknessdata;
		Raw *input = new Raw(l, m, n, inputo);
		/*Filter *f=new Filter();*/
		//input=f->guass3DFilter(input,3);
		RawImage *write = new RawImage();
		ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
		//20140405 delete because of the existance of 
		ls->initialg(*input);
		seedGrowing(input,seedlist);//fisrt time delete the 
		for (int i = 0; i<input->getXsize(); i++)
		{
			for (int j = 0; j<input->getYsize(); j++)
			{
				for (int k = 0; k<input->getZsize(); k++)
				{
					PIXTYPE val = input->get(i, j, k);
					//if (((k-256)*(k-256)+(j-256)*(j-256) )<(230*230))//k<409 && k> 107 && j>156 &&j <390
					//{
					if (val !=1)
					{
						val = 200;  //change to 100 for roc computing *val=0; 
						

					}
					else val = 0; ////change to 0 for roc computing *val=100; 
					input->put(i, j, k, val);
					//}
					//else *val = 0;
				}
			}
		}
		//Seed *seed3 = new Seed(169,367,11);//250,208,12
		seedlist.clear();
		//seedlist.push_back(*seed3);
		//seedGrowing(input,seedlist,0);
		//Seed *seed = new Seed(248, 209, 18);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18:
		//seedlist.push_back(*seed);
		//seedGrowing(input,seedlist,0);
		//ls->initialg(*input);
		seedlist.clear();
		//Seed *seed2 = new Seed(169, 367, 11);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18://250, 208, 12
		seedlist.push_back(*seed);
		Raw *initialdata = initialRegion(seedlist,1, l, m, n);
		*initial=ls->PolypEnergy(*initialdata, *input, 3.0, 0.1, -1, 1.5, 1, iter_outer, pt);//alfa is the ballon pressure,lambda control the curvature
		char *outname1="inner5-8_2polypmethod2.raw";
		char outdir[200]=output;
		strcat(outdir,dirbody);
		strcat(outdir,outname1);
		//for (int i=0;i<input->size();i++)
		//{
		//	if (input->getXYZ(i)!=200)
		//	{
		//		input->putXYZ(i,0);
		//	}
		//}
		for (size_t i = 0; i < initial->getXsize(); i++)
		{
			for (size_t j = 0; j < initial->getYsize(); j++)
			{
				for (size_t k = 0; k < initial->getZsize(); k++)
				{
					if ((i - seed->x)*(i - seed->x) + (j - seed->y)*(j - seed->y) + (k - seed->z)*(k - seed->z) >= 20 * 20)
						initial->put(i, j, k, 0);


				}

			}

		}
		test.writeImageName(*initial,outdir);
	}
	else if (method == 3)//input data: inner data,outer data,no circle data,region growing stop when
	{


		


		
	}
	delete[] indata;

}
float computeVolume(Raw *polyp)
{
	float volume = 0;
	for (size_t i = 0; i < polyp->size(); i++)
	{
		if (polyp->getXYZ(i) != 0)
		{
			volume++;
		}
		

	}
	return volume;
}
float computeArea(Raw *polyp)
{
	ThreeDim_LevelSet *test = new ThreeDim_LevelSet();
	test->initialg(*polyp);
	float area = 0;
	for (size_t i = 0; i < polyp->size(); i++)
	{
		if (polyp->getXYZ(i) != 0)
		{
			area++;
		}


	}
	return area;
}
float computeRadius(Raw * polypg)
{	//min ball contains the set
	/*
	1\compute the gravity center
	2\
	*/
	int x = 0;
	int y = 0;
	int z = 0;
	int n = 0;
	for (size_t i = 0; i < polypg->getXsize(); i++)
	{
		for (size_t j = 0; j < polypg->getYsize(); j++)
		{
			for (size_t k = 0; k < polypg->getZsize(); k++)
			{
				if (polypg->get(i,j,k)!=0)
				{
					n++;
					x += i;
					y += j;
					z += k;





				}

			}

		}//for j..

	}//for i..
	x = x / n;
	y = y / n;
	z = z / n;
	bool flag = true;
	int r = 2;
	while (flag)
	{
		int i = 0;
		int j = 0;
		int k = 0;
		for (size_t i = 0; i < polypg->getXsize() && flag; i++)
		{
			for (size_t j = 0; j < polypg->getYsize() && flag; j++)
			{
				for (size_t k = 0; k < polypg->getZsize() && flag; k++)
				{
					if ((i - x)*(i - x) + (j - y)*(j - y) + (k - z)*(k - z) > r*r)
					{
						flag = false;

					}


				}

			}

		}//for ..i
		if (i==polypg->getXsize())
		{

		}


	}//while








}
void Polyp::computePolyp(Raw *polyp)
{

}
Seed *findplanepoint(Seed *outerPoint,Seed * centerPoint,int height)
{
	
	float gradient = sqrt((outerPoint->x - centerPoint->x) + (outerPoint->y - centerPoint->y) + (outerPoint->z - centerPoint->z));
	float deltax = float(outerPoint->x - centerPoint->x) / gradient;
	float deltay = float(outerPoint->y - centerPoint->y) / gradient;
	float deltaz = float(outerPoint->z - centerPoint->z) / gradient;
	int x = centerPoint->x + deltax*height;
	int y = centerPoint->y+deltay*height;
	int z = centerPoint->z + deltaz*height;
	Seed *planePoint = new Seed(x,y,z);
	return planePoint;
}
void findupplane(Raw *innerandouter,Seed *planePoint,Seed *normaldirection)
{
	vector<Seed> seedlist;
	int x = planePoint->x;
	int y = planePoint->y;
	int z = planePoint->z;
	for (size_t i = 0; i < innerandouter->getXsize(); i++)
	{
		for (size_t j = 0; j < innerandouter->getYsize(); j++)
		{
			for (size_t k = 0; k < innerandouter->getZsize(); k++)
			{
				if (innerandouter->get(i,j,k)!=0)
				{
					if ((i-x)*normaldirection->x+(j-y)*normaldirection->y+(k-z)*normaldirection->z==0)
					{
						innerandouter->put(i, j, k, 200);
					}
				}
			}
		}
	}
	//return innerandouter;
}
void geometrymethod()
{
	int l = 512, m = 512, n = 50;
	vector<Seed>seedlist;
	RawImage *test = new RawImage();
	int offset = 200;
	Raw *outer, *inner;
	PIXTYPE *outerdata = new PIXTYPE[l*m*n + offset];
	PIXTYPE *innerdata = new PIXTYPE[l*m*n + offset];
	outerdata = test->readStreamfloat("D:\\swfdata20140420res\\polyp\\3036P\\inner\\", &l, &m, &n);//outer wall 
	innerdata = test->readStreamfloat("D:\\swfdata20140420res\\polyp\\3036P\\outer\\", &l, &m, &n);//inner wall
	n = 50;
	outer = new Raw(l, m, n, outerdata);
	inner = new Raw(l, m, n, innerdata);
	Raw innerAndOuter = *outer * 2 + *inner;//boundary data,outer is 2,inner is 1
	ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
	ls->initialg(*outer);
	ls->initialg(*inner);
	Seed * middle = new Seed(169, 237, 11);
	seedlist.clear();
	seedlist.push_back(*middle);
	//next step draw a ball
	vector<Seed> outerlist = probabiltyouterpointleastcost(outer, middle, 10);
	Seed *leastcostseed = findLeastCostSeed(outerlist, middle);
	Seed *planePoint = findplanepoint(leastcostseed, middle, 8);//the height is the size of the polyp
	Seed *normal = new Seed(middle->x-planePoint->x,middle->y-planePoint->y,middle->z-planePoint->z);
	findupplane(&innerAndOuter,planePoint,normal);
	vector<Seed> start;
	Seed *seed = new Seed(169,237,11);
	start.push_back(*seed);
	seedGrowing(&innerAndOuter,start);




	
}
void seedGrowing(Raw * src,vector<Seed> seedlist)
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();

	while (!seedqueue.empty())
	{
		
		Seed start=seedqueue.front();
		seedqueue.pop();
		//int size=1;
		//for (int i=start.x-1; i <start.x+1;i++)
		//{
		//	for (int j=start.y-1; j<start.y+1;j++)
		//	{
		//		for (int k=start.z-1; k<start.z+1; k++)
		//		{
		//			if (src->get(i,j,k)==1)
		//			{
		//				Seed *adj=new Seed(i,j,k);
		//				seedqueue.push(*adj);
		//				src->put(i,j,k,200);
		//			}
		//		}
		//	}									
		//}
		int i=start.x;
		int j=start.y;
		int k=start.z;
		//for (int i = 0; i < src->size(); i++)
		//{
		//	src->putXYZ(i,0);

		//}
		if (i>1&&i<src->getXsize()-1&&j>1 && j<src->getYsize()-1 && k > 1&& k<src->getZsize()-1)
		{
			int up=src->get(i-1,j,k)==1;
			int down=src->get(i+1,j,k)==1;
			int left=src->get(i,j-1,k)==1;
			int right=src->get(i,j+1,k)==1;
			int front =src->get(i,j,k-1)==1;
			int back=src->get(i,j,k+1)==1;
			if (src->get(i,j,k)==1)
			{
				int count=0;
				if (up)
				{
					Seed *adj=new Seed(i-1,j,k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj=new Seed(i+1,j,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj=new Seed(i,j-1,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj=new Seed(i,j+1,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj=new Seed(i,j,k-1);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj=new Seed(i,j,k+1);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
					src->put(i,j,k,200);
				//}
				//else src->put(i,j,k,300);
			}
		}


	}//...while


}
void seedGrowing(Raw * src, vector<Seed> seedlist,int threshold)//threshold is the value to be put into the region
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();
	int value = threshold;
	while (!seedqueue.empty())
	{

		Seed start = seedqueue.front();
		seedqueue.pop();

		int i = start.x;
		int j = start.y;
		int k = start.z;

		if (i>1 && i<src->getXsize() - 1 && j>1 && j<src->getYsize() - 1 && k > 1 && k<src->getZsize() - 1)
		{
			int up = src->get(i - 1, j, k) == value;
			int down = src->get(i + 1, j, k) == value;
			int left = src->get(i, j - 1, k) == value;
			int right = src->get(i, j + 1, k) == value;
			int front = src->get(i, j, k - 1) == value;
			int back = src->get(i, j, k + 1) == value;
			if (src->get(i, j, k) == value)
			{
				int count = 0;
				if (up)
				{
					Seed *adj = new Seed(i - 1, j, k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj = new Seed(i + 1, j, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj = new Seed(i, j - 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj = new Seed(i, j + 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj = new Seed(i, j, k - 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj = new Seed(i, j, k + 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
				src->put(i, j , k, 200);
				//}
				//else src->put(i,j,k,300);
			}
		}


	}//...while


}
void seedGrowing(Raw * src,Raw *origin, vector<Seed> seedlist)
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();

	while (!seedqueue.empty())
	{

		Seed start = seedqueue.front();
		seedqueue.pop();
		//int size=1;
		//for (int i=start.x-1; i <start.x+1;i++)
		//{
		//	for (int j=start.y-1; j<start.y+1;j++)
		//	{
		//		for (int k=start.z-1; k<start.z+1; k++)
		//		{
		//			if (src->get(i,j,k)==1)
		//			{
		//				Seed *adj=new Seed(i,j,k);
		//				seedqueue.push(*adj);
		//				src->put(i,j,k,200);
		//			}
		//		}
		//	}
		//}
		int i = start.x;
		int j = start.y;
		int k = start.z;
		int k2 = start.z + 70;
		//for (int i = 0; i < src->size(); i++)
		//{
		//	src->putXYZ(i,0);

		//}
		if (i>1 && i<src->getXsize() - 1 && j>1 && j<src->getYsize() - 1 && k > 1 && k<src->getZsize() - 1)
		{
			int value = 200;
			int originflag = origin->get(i,j,k)>-100;
			int up = src->get(i - 1, j, k) != value && origin->get(i-1, j, k2)>-100;
			//int val = origin->get(i - 1, j, k);
			int down = src->get(i + 1, j, k) != value && origin->get(i + 1, j, k2)>-100;
			int left = src->get(i, j - 1, k) != value && origin->get(i, j - 1, k2)>-100;
			int right = src->get(i, j + 1, k) != value && origin->get(i, j + 1, k2)>-100;
			int front = src->get(i, j, k - 1) != value && origin->get(i, j, k2 - 1)>-100;
			int back = src->get(i, j, k + 1) != value && origin->get(i, j, k2 + 1)>-100;
			if (src->get(i, j, k) != 200)
			{
				int count = 0;
				if (up)
				{
					Seed *adj = new Seed(i - 1, j, k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj = new Seed(i + 1, j, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj = new Seed(i, j - 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj = new Seed(i, j + 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj = new Seed(i, j, k - 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj = new Seed(i, j, k + 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
					src->put(i, j, k, 200);
				//}
				//else src->put(i, j, k, 300);
			}
		}


	}//...while


}
vector<Seed> seedlistdata(const int size)
{
	const int mun = size;
	//data start at 3035P
	int data[][3] =
	{
		386, 267, 187,//3035P sessile
		110, 263, 198,//3035S sessile
		164, 373, 209,//3036P Peduncated
		109, 171, 151,//3036S Peduncated
		191, 289, 157,//3036 P SESSILE
		206, 204, 145,//3036 S SESSILE
		230, 180, 51,//3037P Peduncated
		247, 371, 63,//3037S Peduncated
		356, 286, 201,//3038P Annular
		141, 228, 194,//3038S Annular
		230, 321, 183,//3039P P
		288, 225, 166,//3039s p
		412, 224, 234,//3039P s
		81, 243, 164,//3039S s
		206, 237, 72,//3041P P
		290, 276, 77,//3041S P
		330, 329, 319,//3042P S
		354, 167, 309,//3042S S
		382, 177, 323,	//3042P s
		153, 255, 325,// 3042S 	 s
		234, 155, 318,//3042P s
		156, 283, 327,//3042S s
		94, 262, 362, //3043P s
		393, 228, 354,//3043S s
		337, 219, 294,//3043p S
		137, 286, 288,//3043s S
		300, 350, 107,//3044p S
		195, 196, 102,//3044S s
		211, 372, 175,//3046P P
		300, 177, 147,//3046S p
		397, 204, 218,//3048P p
		63, 231, 177,//3048s S
		142, 195, 224,//3052P P
		414, 255, 216,//3052S P 
		222, 367, 201,//3070p p
		273, 231, 212,//3070s p
		212, 182, 98,//3072P P
		289, 319, 115,//3072S p
		301, 235, 116,//3072P P
		201, 287, 122,//3072S p
		115, 300, 142,// 3072P p
		387, 213, 127,//3072S p
















	};
	vector<Seed> seedlist;
	for (size_t i = 0; i < size; i++)
	{
		seedlist.push_back(Seed(data[i][0], data[i][1], data[i][2]));
	}
	return seedlist;
}
//void ADO()
//{
//	//最好是要把它整理成一个类,要不用起来非常麻烦! 
//	//头文件必须的加下载ADO COM 
//
//		rename("EOF", "adoEOF")//把"EOF"定义改为"adoEOF",避免和文件未eof冲突 
//		//上面头文件的请照写,要不就出现很多错误 
//		//msado15.dll加载后会在我们c++工程项目里自动产生几个头文件,定义了下面 莫名其妙 的宏 类型 等等 
//		inline void sql()
//	{
//		CoInitialize(NULL); //初始化com 
//		_ConnectionPtr m_pConnection = NULL;//数据库连接实例的智能指针 
//		m_pConnection.CreateInstance(__uuidof(Connection));//创建连接实例 
//		try
//		{
//			//参数 1 数据库驱动程序和路径 2 用户名 3 密码 4 照写还不知道是什么 
//			m_pConnection->Open("Provider=Microsoft.Jet.OLEDB.4.0;Data Source=C:\\Database1.mdb",
//				"", "", adModeUnknown);//打开数据库 
//			//也可使用conn->Open("DSN=TSdata","","",adModeUnknown); 
//		}
//		catch (_com_error e)
//		{
//			MessageBox(NULL, "数据库连接失败，确认数据库是否在当前路径下", TEXT("显示窗口"), NULL);
//			return;
//		}
//
//
//		_RecordsetPtr RecordsetPtr; // 创建数据集： 
//
//		try {
//			_bstr_t st = "select * from table1";//com经常用的内部字符串是sql语句表名table1 查所有字段 
//			_variant_t   vRecsAffected; //不知道干什么照写,反正是需要的参数 
//			RecordsetPtr = m_pConnection->Execute(st,             // 利用连接实例和上面2个参数创建数据集： 
//				&vRecsAffected,
//				adOptionUnspecified);//照写 
//		}
//		catch (_com_error *e)
//		{
//			MessageBox(NULL, "创建数据集失败", TEXT("显示窗口"), NULL);
//
//		}
//		while (!RecordsetPtr->GetadoEOF())//数据集的游标如果不是到尾 就执行循环读数据 
//		{
//			_variant_t   vR("ID");//设置要读数据的表的字段名称,这里我哪个数据库里是ID 
//			_variant_t sts = RecordsetPtr->GetCollect(&vR);//取得ID字段的数据 
//
//			_bstr_t st1 = (_bstr_t)sts;//把类型转换一下,要不都不肯显示 
//
//			_variant_t   vR2("字段1");//同上面这个字段是"字段1" 
//			_variant_t sts2 = RecordsetPtr->GetCollect(&vR2);
//			_bstr_t st2 = (_bstr_t)sts2;
//			MessageBox(NULL, (LPCSTR)st1, TEXT("显示ID"), NULL);
//			MessageBox(NULL, (LPCSTR)st2, TEXT("显示字段1"), NULL);
//			RecordsetPtr->MoveNext();//游标移动到下一条记录,准备下次循环 
//		}
//		CoUninitialize();//释放com 
//	}
//}