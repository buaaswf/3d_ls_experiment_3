
#include "vol_math_LevelSet.h"
//#include "statistics.h"
#include "test.h"
#include "ProcessDirty.h"
#include "DeleteMiddle.h"
#include"Polyp.h"
#include"GoodData.h"
#include"Filter.h"
#ifdef _WIN32
#include <Windows.h>
#include <strsafe.h>
#else
#include <dirent.h>
#endif
#define output "D:\\swfdata20140420res\\polyp\\" 
//#define  input2 "D:\\data\\clean\\polypseginputdata\\origin\\"
//#define output "K:\\20140404\\" 
//#define input1  "L:\\sdfdata2\\edt\\20140409edt\\"	//thickness uint8	//edt		//float
////swf 20140409 delete for float2char
//#define input2  "K:\\sdf\\volume\\clean\\clean\\ep\\20140410\\"//short
////#define input2 "K:\\20140404\\inner\\"
//#define input3  "K:\\skeleton\\"  //skeleton uint8 //unsigned char
//#define input1  "F:\\data\\skeleton-edt\\"				//float
//#define input2	"E:\\volume\\skeletono\\"		//short
//#define input3 "F:\\data\\skeleton-s\\"   //unsigned char
//#define input1 "L:\\sdfdata2\\edt\\20140414edt\\"

//#define input1 "D:\\swfdata20140420res\\auto\\outer\\"//GAC levelset data
//# define input1 "L:\\sdfdata2\\inner\\"
//#define input2 "E:\\volume\\segmention\\rate\\"
//# define input2 "F:\\data\\dirty\\bone\\"K:\sdf\volume\clean\clean\ep\clean
//#define input2 "K:\\sdf\\volume\\clean\\clean\\ep\\clean\\"
// roc计算GAC的FP.TP值用的
// #define input1 "D:\\swfdata20140420res\\roc3041\\"
//#define input2 "D:\\segdata\\people\\roc3041\\" //people data roc 3041
//#define input3 "K:\\sdf\\volume\\clean\\clean\\ep\\clean\\roc3041\\"
//#define  input3 "L:\\sdfdata2\\edt\\20140414skeleton\\"
//#define input1 "D:\\swfdata20140420res\\auto\\edt\\"
//#define input2 "K:\\sdf\\volume\\clean\\clean\\ep\\clean\\"
//#define input3 "D:\\swfdata20140420res\\auto\\skeleton\\"
/*
35-44divide region dir
#define input1 "L:\\sdfdata2\\edt\\35-44edt\\"
#define input2 "K:\\sdf\\volume\\clean\\clean\\ep\\clean\\"
#define input3 "K:\\skeleton\\35-44\\"
*/
//people divide region dir
//#define input1 "D:\\segdata\\people\\thickness\\edt\\"
#define input1 "L:\\sdfdata2\\inner\\"
//#define input2 "K:\\sdf\\volume\\clean\\clean\\ep\\clean\\"
#define input3 "D:\\segdata\\people\\skeleton\\"
// people /people divide region dir
#define  input2 "D:\\data\\clean\\polypseginputdata\\origin\\"
using namespace cimg_library;
using namespace std;
//////////////////////////////////////////////////////////////////////////
//获取指定目录下所有文件的文件名，不包括文件夹，在GetFileFromDir中使用
//strDir: 输入，目录路径
//vFileDirList： 输出，文件路径列表
//返回：空
//////////////////////////////////////////////////////////////////////////
void GetFileNameFromDirv2(string strDir, vector<string>& vFileDirList)
{
#ifdef _WIN32
	WIN32_FIND_DATAA ffd;
	LARGE_INTEGER filesize;
	string szDir;
	//size_t length_of_arg;
	HANDLE hFind = INVALID_HANDLE_VALUE;
	DWORD dwError = 0;

	szDir = strDir + "\\*";
	hFind = FindFirstFileA(szDir.c_str(), &ffd);

	if (INVALID_HANDLE_VALUE == hFind)
	{
		cout << "get file name error" << endl;
		return;
	}
	do
	{
		if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
		{
			string filename = ffd.cFileName;//(const char*)
			string filedir = strDir + "\\" + filename;
			vFileDirList.push_back(filedir);
		}
	} while (FindNextFileA(hFind, &ffd) != 0);

	dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES)
	{
		cout << "FindFirstFile error" << endl;
		return;
	}
	FindClose(hFind);
#else
	DIR *dir;
	struct dirent *ptr;
	dir = opendir(strDir.c_str());
	while ((ptr = readdir(dir)) != NULL)
	{
		string path = strDir + string("/") + string(ptr->d_name);
		vFileDirList.push_back(path);
	}
	closedir(dir);
	sort(vFileDirList.begin(), vFileDirList.end());
#endif
}

//#include <iostream>
//#include <crtdbg.h> 
//#include "CImg.h" 
//#include "ThreeDim_LevelSet.h"
//#include "Filter.h"
//#include "WipeNioisePde.h"
//#include "string.h"
////================
//#include <fstream>
//#include <vector>
//#include <string>
//#ifdef _WIN32
//#include <Windows.h>
//#include <strsafe.h>
//#else
//#include <dirent.h>
//#endif
//#define output "D:\\sdfdata\\" 
//#define input1  "D:\\sdfdata\\pvaluethickness\\" //"K:\\sdf\\volume\\clean\\clean\\ep\\""E:\\volume\\cleantestdata\\test2\\"
//#define input2  "K:\\sdf\\volume\\clean\\clean\\ep\\test3\\" 
//using namespace cimg_library;
//using namespace std;


void testcolon(int argc,string dir)
{
	
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=10;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	for (int i = 0; i < l*m*n; i++)
	{
		inputo[i]=(float) indata[i];		
	}

	Raw *input=new Raw(l,m,n,inputo);
	Filter *f=new Filter();
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//20140405 delete because of the existance of 
	ls->initialg(*input);
	for (int i=0; i<input->getXsize(); i++)
	{
		for (int j=0; j<input->getYsize(); j++)
		{
			for (int k=0; k<input->getZsize(); k++)
			{
				if (input->get(i,j,k)>=1)
				{
					initial->put(i,j,k,-2);
				}
				else initial->put(i,j,k,2);

			}
		}

	}
	*initial=ls->minimal_surface(*initial,*input,5.0,0.1,-3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//Raw temp(*initial);
	ls->outerwallauto(*initial,*input,10,0.1,-6,1.5,1,10,pt);
	//*initial -=temp;
	char *outname2="outer5-8_2_20140405.raw";
	char outdir2[200]=output;
	strcat(outdir2,dirbody);
	strcat(outdir2,outname2);
	test.writeImageName(*initial,outdir2);
//	evaluate(dir,l,m,n);
}

void testsesmic()
{
	char *pt="single_well";
	int l = 201,m = 201,n = 851, l1=0,l2=0,iter_outer = 10;
	RawImage test;
	char dirbody[100];
	unsigned char * indata=new unsigned char[l*m*n];
	test.readImage(indata,"K:\\sdf\\geo\\Probe_fault_Amp.probe .raw",l*m*n);//F:\\PA1\\ST1\\SE1\\  //K:\\sdf\\MRI
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	short min = 1000,max = -100;
	for (int i = 0; i < l*m*n; i++)
	{
		//change the big --little
		//float * p= (float *)(indata+i);
		//unsigned char * bp= (unsigned char *)p;
		//std:swap(bp[0],bp[3]);
		//std::swap(bp[1],bp[2]);
		min < indata[i] ? min=min:min=indata[i];
		max > indata[i] ? max=max:max=indata[i];
		//cal the max and min data
	/*	if ( indata[i] >= 864 && indata[i] <= 1063 )
		{
			inputo[i] = 100;
		} 
		else
		{
			inputo[i] = (short )0;
		}
		*/
		inputo[i]=(float) indata[i];		
	}

	cout <<min << max <<endl;

	Raw *input=new Raw(l,m,n,inputo);

	//Filter *f=new Filter();
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//ls->initialg(*input);
	for (int i=0; i<input->getXsize(); i++)
	{
		for (int j=0; j<input->getYsize(); j++)
		{
			for (int k=0; k<input->getZsize(); k++)
			{
				//if (input->get(i,j,k) >= 1)
				//{
				//	initial->put(i,j,k,-2);
				//}
				//else 
					//if ((i >= 172 && i <= 352 && j >= 164 && j <= 376 && z>19 && z <))
				if ((i >= 196 && i <= 220 && j >= 202 && j <= 267 && k > 40 && k < 50))
				{
					initial->put(i, j, k, -2);
				} 
				else
				{
					initial->put(i, j, k, 2);
				}


			}
		}

	}
	*initial=ls->minimal_surface(*initial,*input,5.0,0.1,-3,1.5,1,iter_outer,pt);
	//if you available this, don,t
	//forget to change the next line to initial
	test.writeMRI(*initial,"K:\\sdf\\geo\\data.raw");//F:\\PA1\\ST1\\SE1

}
//void testhistgram()
//{
//	//HUandThickness();
//	directdivideregion();
//}
void rate(string dir)
{
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=10;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	delete [] indata;
	
}
void makeGoodData()
{
	for (int i=0;i<10;i++)
	{
		//double d=(double)i;
		Raw *data=myColondata4(512,512,50,0.9+i*0.01);
		RawImage *outdata=new RawImage();
		char file_no[4];
		int filen = i;
		string strDir("D:\\goodata0.9");
		itoa(filen, file_no, 10);//把数字存储为char的数组
		strDir += file_no;//string是标准库类型，可以直接与char的数组进行+号连接
		strDir += ".raw";
		const char *cstr = strDir.c_str();
		outdata->writenormal(*data,cstr);
		delete data;
	}
}


//void polypseg()
//{
//	string dir2(input2);//D:\swfdata20140420res\polyp\regiongrowingfillnull
//	string dirthickness(inputt);
//	vector<string> files2;
//	vector<string> filesthickness;
//	GetFileNameFromDir(dir2, files2);
//	GetFileNameFromDir(dirthickness, filesthickness);
//
//	vector<string> ::iterator thicknessiter = filesthickness.begin();
//	vector<string>::iterator iterFile2;
//	int i = 0;
//	for (iterFile2 = files2.begin(); iterFile2 != files2.end(); iterFile2++)
//	{
//
//		i++;
//		iterFile2->assign(iterFile2->substr(dir2.size() + 1));
//		thicknessiter->assign(thicknessiter->substr(dirthickness.size() + 1));
//		cout << *iterFile2 << endl;
//		Polyp *test = new Polyp();
//		test->polypDetect(*iterFile2, *thicknessiter, 2,i);
//		thicknessiter++;
//	}
//
//
//}
int main(int argc,char **argv)
{
	string dir2(input2);//D:\swfdata20140420res\polyp\regiongrowingfillnull
	string dirthickness(inputt);
	string dirseg("D:\\data\\segmention\\polypseg\\");
	vector<string> files2;
	vector<string> filesthickness;
	vector<string> fileseg;
	GetFileNameFromDirv2(dir2,files2);
	GetFileNameFromDirv2(dirthickness,filesthickness);
	GetFileNameFromDirv2(dirseg, fileseg);
	//seedlistdata();
	int cur = 1;
	vector<string> ::iterator thicknessiter = filesthickness.begin() + cur-1;
	vector<string> ::iterator segdiriter = fileseg.begin() + cur-1;
	vector<string>::iterator iterFile2;
	int i = 0;
	for ( iterFile2 = files2.begin()+cur-1; iterFile2 != files2.end(); iterFile2++ )
	{

		
		iterFile2->assign(iterFile2->substr(dir2.size()+1));
		thicknessiter->assign(thicknessiter->substr(dirthickness.size()+1));
		segdiriter->assign(segdiriter->substr(dirseg.size() + 1));
		cout<<*iterFile2 <<endl;
		Polyp *test=new Polyp();

		test->polypDetect(*iterFile2, *thicknessiter, *segdiriter, 2, i+cur-1);
		i++;
		segdiriter++;
		thicknessiter++;
		//ddcircle(*iterFile);
		//testcolon(argc,*iterFile2);
		//float2uchar(512,512,700,*iterFile2);
		//testsesmic();
		//thincknessstdv2(*iterFile2);
		//roc(*iterFile2,);
		//rate(*iterFile2);

	}
	//testhistgram();
	//deleteMiddle(argc,"");
	//test();//delete dirty success
	//cout<<endl;
	//roc3();
	//threshold();
	//rocway2();
	//testcolontest();
	//
	//testsesmic();

	
	system("pause");
	return 0;

}
