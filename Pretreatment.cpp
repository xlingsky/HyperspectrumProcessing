// Pretreatment.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"
#include "Pretreatment.h"
#include <math.h>
#include <string.h>
#include <malloc.h>

#ifndef _NO_IMAGEIO
#include "package_gdal/imagebase.h"
#include "DPGRID/SpWallisFlt.hpp"

class CWuWlsImage : public CSpWlsImage
{
public:
	CWuWlsImage(){}
	~CWuWlsImage(){
		m_img.Close();
	}
	bool Open(const char* lpstrPathName){
		return m_img.Open(lpstrPathName);
	}
	bool Create(const char* lpstrPath, int nCols, int nRows){
		return m_img.Create(lpstrPath, nCols, nRows, 1);
	}

	virtual int GetRows() { return m_img.GetRows(); }
	virtual int GetCols() { return m_img.GetCols(); }
	virtual int GetPixelBytes() { return 1; }
	virtual BOOL Read(BYTE* pBuf, int rowIdx){
		return m_img.ReadBand8(pBuf,0,0,rowIdx,m_img.GetCols(),1)?TRUE:FALSE;
	}
	virtual BOOL Write(BYTE *pBuf, int rowIdx){
		return m_img.Write(pBuf, 0, rowIdx, m_img.GetCols(), 1)?TRUE:FALSE;
	}
	virtual void GetFilePath(char* strPath){
		strcpy(strPath, m_img.GetImagePath());
	}
	CImageBase* GetImageBase()	{ return &m_img; }
	void		CopyGeoInformation(CWuWlsImage& img){
		m_img.CopyGeoInformation(*(img.GetImageBase()));
	}
private:
	CImageBase m_img;
};

inline void GetWallisParFile(char* strPath){
	char* pS = strrchr(strPath, '.');	if (pS) strcpy(pS, ".wallis"); else strcat(strPath, ".wallis");
}

class CWuWlsFile : public CSpWlsFile
{
public:
	BOOL WallisFlt(const char* lpstrImagePath, const char* lpstrWallisPath){
		CWuWlsImage imgR, imgW;
		if (!imgR.Open(lpstrImagePath)) return false;
		int nRows = imgR.GetRows();
		int nCols = imgR.GetCols();
		if (!imgW.Create(lpstrWallisPath, nCols, nRows)) return false;
		imgW.CopyGeoInformation(imgR);
		return WallisFlt(&imgR, &imgW) ;
	}
	BOOL WallisFlt(CWuWlsImage *pImgFileS, CWuWlsImage *pImgFileD){
		int colsS = pImgFileS->GetCols(), rowsS = pImgFileS->GetRows(), pxlBytes = pImgFileS->GetPixelBytes();
		
		char strPath[512];		pImgFileS->GetFilePath(strPath);
		{
			GetWallisParFile(strPath);
			if (!Load4File(strPath,TRUE)){
				if (!CalWlsPar(pImgFileS)) return FALSE;
				Save2File(strPath);
			}
		}

		PrintMsg("Apply WallisFilter... ");		ProgBegin(rowsS / 512+1);	int cancel = 0;

		int lineSize = colsS*pxlBytes, vv;
		BYTE *pRow = new BYTE[lineSize], *pD, *pS;
		for (int r = 0; r<rowsS; r++){
			pImgFileS->Read(pRow, r);
			Bnd2RGB(pRow, colsS, pxlBytes);
			if (pxlBytes>1){ for (pD = pRow, pS = pRow, vv = 0; vv < colsS; vv++, pD++, pS += 3){ *pD = BYTE((UINT(*pS) + *(pS + 1) + *(pS + 2)) / 3); } }

			WlsFlt(pRow, r, colsS);
			pImgFileD->Write(pRow, r);

			if ((r % 512) == 0) ProgStep(cancel);
		}
		delete pRow;

		ProgEnd();
		return TRUE;
	}

protected:
	virtual void ProgBegin(int range)       { printf("%5d/%5d", range, 0);	m_step = 0; m_range = range; };
	virtual void ProgStep(int& cancel)      { m_step++;	printf("\b\b\b\b\b%5d", m_step); };
	virtual void ProgEnd()                  { printf("\b\b\b\b\b%5d\n", m_range); };
	virtual void PrintMsg(LPCSTR lpstrMsg) { printf("%s",lpstrMsg); };
private:
	int m_step;
	int m_range;
};
bool Pretreatment::Wallisfilter(const char* lpstrImagePath, const char* lpstrWallisPath){
	CWuWlsFile wlsFile;
	return wlsFile.WallisFlt(lpstrImagePath, lpstrWallisPath) ? true : false;
	//	return ::Wallisfilter(lpstrImagePath,lpstrWallisPath, WALLIS_meanValue, WALLIS_sigmaValue, WALLIS_fC_Value, WALLIS_fB_Value);
}
bool Pretreatment::CheckWallisParConfig(const char* lpstrImagePath){
	char strPath[512];	strcpy(strPath, lpstrImagePath);
	GetWallisParFile(strPath);
	CSpWlsFile wlsFile;
	return wlsFile.Load4File(strPath, TRUE) ? true : false;
}
#endif

#define WALLIS_gridWIndowSize 16
#define WALLIS_filterWindowSize 31	
#define WALLIS_fC_Value 0.8f 
#define WALLIS_fB_Value  0.9f		//0.85f 
#define WALLIS_meanValue  137.0f	//127.0f 
#define WALLIS_sigmaValue  190.0f	//60.0f  

#define PI 3.1415926

inline bool	CalWallisParameter(float *r0, float *r1,
	BYTE *data, int winPtSum,
	float meanValue, float sigmaValue,
	float C_Value, float B_Value)
{
	int	i;
	float mean, sigma;
	BYTE	*pData;

	mean = sigma = 0.0;
	pData = data;
	for (i = 0; i < winPtSum; i++)
	{
		mean += *pData;
		sigma += *pData * *pData;
		pData++;
	}

	mean = mean / winPtSum;
	sigma = sigma / winPtSum - mean*mean;
	if (sigma < 0) sigma = 0;
	sigma = float(sqrt(sigma));

	if (sigma == 0.0)
	{
		*r1 = 1.0;
		*r0 = B_Value*meanValue + (1.0f - B_Value - *r1)*mean;
	}
	else
	{
		*r1 = C_Value*sigmaValue / (C_Value*sigma + (1.0f - C_Value)*sigmaValue);
		*r0 = B_Value*meanValue + (1.0f - B_Value - *r1)*mean;
	}

	return 1;
}
// 
// inline bool	CalWallisParameter(float *r0, float *r1,
// 	CImageBase& img, int stCol, int stRow, int nCols, int nRows,
// 	float meanValue, float sigmaValue,
// 	float C_Value, float B_Value)
// {
// 	int	i,j;
// 	float mean, sigma;
// 	BYTE data;
// 	int winPtSum = nCols*nRows;
// 
// 	mean = sigma = 0.0;
// 	
// 	for (i = 0; i < nRows; i++)
// 	{
// 		for (j = 0; j < nCols; j++){
// 			data = img.GetBandVal8(stCol+j,stRow+i,0);
// 			mean += data;
// 			sigma += data*data;
// 		}
// 		
// 	}
// 
// 	mean = mean / winPtSum;
// 	sigma = sigma / winPtSum - mean*mean;
// 	if (sigma < 0) sigma = 0;
// 	sigma = float(sqrt(sigma));
// 
// 	if (sigma == 0.0)
// 	{
// 		*r1 = 1.0;
// 		*r0 = B_Value*meanValue + (1.0f - B_Value - *r1)*mean;
// 	}
// 	else
// 	{
// 		*r1 = C_Value*sigmaValue / (C_Value*sigma + (1.0f - C_Value)*sigmaValue);
// 		*r0 = B_Value*meanValue + (1.0f - B_Value - *r1)*mean;
// 	}
// 
// 	return 1;
// }


inline float	 InterplotWallisParameter(int x, int y, float *grid, int gridCol, int gridRow)
{
	int      grid_r, grid_c;

	grid_r = y / WALLIS_gridWIndowSize;
	grid_c = x / WALLIS_gridWIndowSize;

	if (grid_r < 0 || grid_c < 0 || grid_r >= gridRow - 1 || grid_c >= gridCol - 1)
	{
		if (grid_r <= 0)         grid_r = 0;
		if (grid_c <= 0)         grid_c = 0;

		if (grid_r >= gridRow - 1) grid_r = gridRow - 1;
		if (grid_c >= gridCol - 1) grid_c = gridCol - 1;

		x = grid_c*WALLIS_gridWIndowSize;
		y = grid_r*WALLIS_gridWIndowSize;
	}

	float *pBuf, ratioX, ratioX1, ratioY, ratioY1;
	ratioX = float(x - grid_c*WALLIS_gridWIndowSize) / WALLIS_gridWIndowSize;
	ratioX1 = 1 - ratioX;
	ratioY = float(y - grid_r*WALLIS_gridWIndowSize) / WALLIS_gridWIndowSize;
	ratioY1 = 1 - ratioY;

	pBuf = grid + grid_r*gridCol + grid_c;

	float	XZ1 = pBuf[0] * ratioX1 + pBuf[1] * ratioX;
	pBuf += gridCol;
	float	XZ2 = pBuf[0] * ratioX1 + pBuf[1] * ratioX;
	float	YZ = XZ1*ratioY1 + XZ2*ratioY;

	return YZ;
}

inline void pick_buffer_block(BYTE *pSrc, int nCols, int nRows,
	int stCol, int stRow, BYTE* pDst, int winCol, int winRow)
{
	int i;
	if (stCol < 0) stCol = 0;	if (stRow < 0) stRow = 0;
	BYTE* pBuf = pSrc + stRow*nCols;
	for (i = 0; i < winRow; i++)
	{
		memcpy(pDst, pBuf + stCol, winCol);
		pBuf += nCols;	pDst += winCol;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
//CPretreatment

// ////////////////////////////////////////////////////////////////////////////////////////
// /*
// 			函数名称：TransToGray
// 			函数功能:将RGB影像转化为灰度影像
// 			参数说明：
// 			m_xsize    影像的宽
// 			m_ysize    影像的高
// 			ph_after   输出的灰度影像
// 			p1    原始影像的红波段
// 			p2    原始影像的绿波段
// 			p3    原始影像的蓝波段			
// */
// ////////////////////////////////////////////////////////////////////////////////////////
// void CPretreatment::TransToGray(int m_xsize,int m_ysize,WORD *ph_after,
// 							    WORD *p1,WORD *p2,WORD *p3) //把24位彩色图像转换成8位的灰度图像
// {
// 	int i,j;
// 	
// 	//对源图像进行操作
// 	// 每行 
// 	for(i = 0; i < m_ysize; i++)
// 	{
// 		//每列
// 		for(j = 0;j < m_xsize ;j++)
// 		{
// 			ph_after[i * m_xsize +j] = (WORD)(0.299 * p1[i * m_xsize +j]+
// 										   	  0.587 * p2[i * m_xsize +j]+
// 											  0.114 * p3[i * m_xsize +j]+0.5);
// 		}
// 	}	
// }

//////////////////////////////////////////////////////////////////////////wallis变换
/*
	//wallis 滤波:
	//参考：Wallis 滤波在影像匹配中的应用. 武汉测绘科技大学学报,1999 ,24 (1) :24～27
	float fC_Value =0.8;
	float fB_Value =0.8;
	Wallisfilter(imHeight,imWidth,			//图像高宽
	127.0f,					//目标图像灰度的平均值---一般为127
	60,						//目标图像的标准差---一般为40～65
	fC_Value,				//fC_Value ∈[0 , 1 ]	为影像方差的扩展常数， 取值范围为0. 75～1 
	fB_Value,				//fB_Value ∈[0 , 1 ]  为影像的亮度系数， 取值范围为0. 5～1
	pData);					//图像灰度数据

*/
//
void Pretreatment::Wallisfilter(BYTE* data, int nCols, int nRows)
{
	Wallisfilter(data,nCols,nRows,WALLIS_meanValue,WALLIS_sigmaValue,WALLIS_fC_Value,WALLIS_fB_Value); 
}

void Pretreatment::Wallisfilter(BYTE* data, int nCols, int nRows,
				  float meanValue,float sigmaValue,
				  float C_Value,float B_Value ) 
{      
	int   i,j;	
	int	  x,y,br,bc,gridRow,gridCol;
	float rc,rmean,rsigma,*gridR0,*gridR1;
	BYTE *block;

	gridRow = (nRows-WALLIS_filterWindowSize)/WALLIS_gridWIndowSize;
	gridCol = (nCols-WALLIS_filterWindowSize)/WALLIS_gridWIndowSize;
	
	block	= (BYTE *)calloc(WALLIS_filterWindowSize*WALLIS_filterWindowSize, sizeof(BYTE));
	gridR0 = (float *)calloc( (gridRow+1)*(gridCol+1), sizeof(float) );
	gridR1 = (float *)calloc( (gridRow+1)*(gridCol+1), sizeof(float) );
	
	rmean = rsigma = 0.0;
	BYTE *pPos = data;   float *pR0 = gridR0, *pR1 = gridR1;
	BYTE *pCenter = block + WALLIS_filterWindowSize*WALLIS_filterWindowSize / 2;
	br = 0;
	for( i=0; i<gridRow; i++ )
	{		
		bc = 0;
		for (j=0; j<gridCol; j++) 
		{
			pick_buffer_block(pPos, nCols, 0, bc, 0, block, WALLIS_filterWindowSize, WALLIS_filterWindowSize);
			CalWallisParameter( pR0,pR1,block, WALLIS_filterWindowSize*WALLIS_filterWindowSize,meanValue, sigmaValue, C_Value, B_Value );
			
			rc = *pCenter**pR1 + *pR0;
			rmean += rc;
			rsigma += rc*rc;

			bc += WALLIS_gridWIndowSize;
			pR0++;	pR1++;
		}
		br		+= WALLIS_gridWIndowSize;
		pPos += WALLIS_gridWIndowSize*nCols;
	}
	free(block);
	
	rc = (float)gridRow*gridCol;
	rmean = rmean/rc;
	rsigma = rsigma/rc-rmean*rmean;
	if (rsigma<0) rsigma = 0;
	rsigma = float( sqrt(rsigma) );
	
	float r0, r1;	int gf;
	pPos = data; 	
	y =  - WALLIS_filterWindowSize / 2;
	for( i=0; i<nRows; i++ ) 
	{
		x =  - WALLIS_filterWindowSize / 2;
		for (j = 0; j<nCols; j++, pPos++)
		{
			
			r0 = InterplotWallisParameter(x, y, gridR0, gridCol, gridRow);
			r1 = InterplotWallisParameter(x, y, gridR1, gridCol, gridRow);
			
			gf = *pPos;
			if (gf<=3 || gf>=252) continue;			
			gf = int((gf*r1+r0-rmean)*sigmaValue/rsigma+127);
			if (gf >= 256) gf = 255; else if (gf<0) gf = 0;
			*pPos = (BYTE)gf;
			x++;
		}		
		y++;
	}
	free(gridR0);	 free(gridR1);
	return  ;
}

// #define MAX_BUFFER_SIZE		(1024*1024*128)
// bool Wallisfilter(const char* lpstrImagePath,const char* lpstrWallisPath,
// 	float meanValue, float sigmaValue,
// 	float C_Value, float B_Value)
// {
// 	CImageBase img,imgW;
// 	if (!img.Open(lpstrImagePath)) return false;
// 
// 	int nRows = img.GetRows();
// 	int nCols = img.GetCols();
// 
// 	if (!imgW.Create(lpstrWallisPath, nCols, nRows, 1)) return false;
// 
// 	printf("----Wallis[%s->%s]----\n", lpstrImagePath, lpstrWallisPath);
// 
//  	BYTE* data = new BYTE[nCols+8];
// // 	int nRows_Read = nRows;		int nRows_Read_start = 0;
// // 	if (nRows*nCols<MAX_BUFFER_SIZE){
// // 		data = new BYTE[nRows*nCols];
// // 		img.Read8(data, 0, 0, nCols, nRows);
// // 	}
// // 	else
// // 	{
// // 		nRows_Read = MAX_BUFFER_SIZE / nCols;
// // 		nRows_Read = nRows_Read / WALLIS_gridWIndowSize *WALLIS_gridWIndowSize;
// // 		data = new BYTE[nCols*nRows_Read];
// // 		img.Read8(data, 0, 0, nCols, nRows_Read);
// // 		nRows_Read_start = nRows_Read;
// // 	}
// 
// 	int   i, j;
// 	int	  x, y, br, bc, gridRow, gridCol;
// 	float rc, rmean, rsigma, *gridR0, *gridR1;
// 	BYTE *block;
// 
// 	gridRow = (nRows - WALLIS_filterWindowSize) / WALLIS_gridWIndowSize;
// 	gridCol = (nCols - WALLIS_filterWindowSize) / WALLIS_gridWIndowSize;
// 
// 	block = (BYTE *)calloc(WALLIS_filterWindowSize*WALLIS_filterWindowSize, sizeof(BYTE));
// 	gridR0 = (float *)calloc((gridRow + 1)*(gridCol + 1), sizeof(float));
// 	gridR1 = (float *)calloc((gridRow + 1)*(gridCol + 1), sizeof(float));
// 
// 	printf("bulid %d*%d r0 r1 grid...\n",gridCol,gridRow );
// 
// 	rmean = rsigma = 0.0;
// //	BYTE *pPos = data;   
// 	float *pR0 = gridR0, *pR1 = gridR1;
// 	BYTE *pCenter = block + WALLIS_filterWindowSize*WALLIS_filterWindowSize / 2;
// 	br = 0;
// 	for (i = 0; i < gridRow; i++)
// 	{
// 		bc = 0;
// 		for (j = 0; j < gridCol; j++)
// 		{
// //			pick_buffer_block(pPos, nCols, 0, bc, 0, block, WALLIS_filterWindowSize, WALLIS_filterWindowSize);
// 			CalWallisParameter(pR0, pR1, img, bc, br, WALLIS_filterWindowSize,WALLIS_filterWindowSize,
// 				 meanValue, sigmaValue, C_Value, B_Value);
// //			CalWallisParameter(pR0, pR1, block, WALLIS_filterWindowSize*WALLIS_filterWindowSize, meanValue, sigmaValue, C_Value, B_Value);
// 
// 			rc = *pCenter**pR1 + *pR0;
// 			rmean += rc;
// 			rsigma += rc*rc;
// 
// 			bc += WALLIS_gridWIndowSize;
// 			pR0++;	pR1++;
// 		}
// 		br += WALLIS_gridWIndowSize;
// // 		if (br>=nRows_Read){
// // 			img.Read8(block, bc, br, WALLIS_filterWindowSize, WALLIS_filterWindowSize);
// // 		}
// // 		pPos += WALLIS_gridWIndowSize*nCols;
// 	}
// 	free(block);
// 
// 	rc = gridRow*gridCol;
// 	rmean = rmean / rc;
// 	rsigma = rsigma / rc - rmean*rmean;
// 	if (rsigma < 0) rsigma = 0;
// 	rsigma = float(sqrt(rsigma));
// 
// 	float r0, r1;	int gf;
// 	printf("convert pixel value...\n");
// //	pPos = data;
// 	y = -WALLIS_filterWindowSize / 2;
// 	for (i = 0; i < nRows; i++)
// 	{
// 		x = -WALLIS_filterWindowSize / 2;
// 		for (j = 0; j < nCols; j++)//, pPos++
// 		{
// 
// 			r0 = InterplotWallisParameter(x, y, gridR0, gridCol, gridRow);
// 			r1 = InterplotWallisParameter(x, y, gridR1, gridCol, gridRow);
// 
// //			gf = *pPos;
// 			gf = img.GetBandVal8(j, i, 0);
// 			if (gf <= 3 || gf >= 252) continue;
// 			gf = int((gf*r1 + r0 - rmean)*sigmaValue / rsigma + 127);
// 			if (gf >= 256) gf = 255; else if (gf < 0) gf = 0;
// //			*pPos = (BYTE)gf;
// 			data[j] = (BYTE)gf;
// 			x++;
// 		}
// 		y++;
// 		imgW.Write(data, 0, i, nCols, 1);
// 	}
// 	free(gridR0);	 free(gridR1);
// 	delete data;
// 	printf("----wallis end----\n");
// 	return true;
// }
// 

inline double* MakeGauss(double sigma, int *pnWindowSize)
{
	// 循环控制变量
	int i;

	// 数组的中心点
	int nCenter;

	// 数组的某一点到中心点的距离
	double  dDis;

	// 中间变量
	double  dValue;
	double  dSum;
	dSum = 0;

	// 数组长度，根据概率论的知识，选取[-3*sigma, 3*sigma]以内的数据。
	// 这些数据会覆盖绝大部分的滤波系数
	nCenter = int(ceil(2 * sigma));
	*pnWindowSize = 1 + 2 * nCenter;

	// 分配内存
	//*pdKernel = new double[*pnWindowSize] ;
	double* pdKernel = (double*)malloc(*pnWindowSize*sizeof(double));/////////////////////
	for (i = 0; i< (*pnWindowSize); i++)
	{
		dDis = (double)(i - nCenter);
		dValue = exp(-0.5*dDis*dDis / (sigma*sigma)) / (sqrt(2 * PI) * sigma);
		pdKernel[i] = dValue;
		dSum += dValue;
	}

	// 归一化
	for (i = 0; i<(*pnWindowSize); i++)
	{
		pdKernel[i] /= dSum;
	}
	return pdKernel;
}

void Pretreatment::GaussianSmooth(unsigned char* data, int nCols, int nRows, double sigma)
{
	// 循环控制变量
	int y;
	int x;

	int i;

	// 高斯滤波器的数组长度

	int nWindowSize;

	//  窗口长度的1/2
	int	nHalfLen;

	// 一维高斯数据滤波器
	double *pdKernel;

	// 高斯系数与图象数据的点乘
	double  dDotMul;

	// 高斯滤波系数的总和
	double  dWeightSum;

	// 中间变量
	double * pdTmp;

	// 分配内存
	pdTmp = (double*)malloc(nCols*nRows*sizeof(double));

	// 产生一维高斯数据滤波器
	pdKernel = MakeGauss(sigma, &nWindowSize);

	// MakeGauss返回窗口的长度，利用此变量计算窗口的半长
	nHalfLen = nWindowSize / 2;

	// x方向进行滤波
	unsigned char* pData = data;
	double* pTmp = pdTmp;
	for (y = 0; y<nRows; y++)
	{
		for (x = 0; x<nCols; x++, pData++, pTmp++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for (i = (-nHalfLen); i <= nHalfLen; i++)
			{
				// 判断是否在图象内部
				if ((i + x) >= 0 && (i + x) < nCols)
				{
					dDotMul += (double)(*(pData + i )) * (*(pdKernel + nHalfLen + i));////  ?
					dWeightSum += *(pdKernel + nHalfLen + i);
				}
			}
			*pTmp = dDotMul / dWeightSum;
		}
	}

	// y方向进行滤波
	for (x = 0; x<nCols; x++)
	{
		pData = data+x;	 pTmp = pdTmp+x;
		for (y = 0; y < nRows; y++, pData += nCols, pTmp += nCols)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for (i = (-nHalfLen); i <= nHalfLen; i++)
			{
				// 判断是否在图象内部
				if ((i + y) >= 0 && (i + y) < nRows)
				{
					dDotMul += (double)(*(pTmp +i*nCols)) * (*(pdKernel + nHalfLen + i));
					dWeightSum += *(pdKernel + nHalfLen + i);
				}
			}
			*(pData) = (int)(dDotMul / dWeightSum);
		}
	}

	// 释放内存
	free(pdKernel);
	free(pdTmp);
}

bool Pretreatment::AffineResample(const BYTE* data, int nCols, int nRows, BYTE* resampleData, int nReCols, int nReRows, int stBlockCol, int stBlockRow, int nBlockCols, int nBlockRows, float aop6[6],RESAMPLEMODE mode )
{
	int r, c;
	int nx, ny;

	float x1, y1, x2, y2;
//	float dxx = aop6[0], dxy = aop6[1], dyx = aop6[3], dyy = aop6[4];

	BYTE* pBuf = resampleData + stBlockRow*nReCols + stBlockCol;
//	x1 = (stBlockCol)*aop6[0] + (stBlockRow )*aop6[1] + aop6[2];
//	y1 = (stBlockCol)*aop6[3] + (stBlockRow )*aop6[4] + aop6[5];
	float offsetX = 0, offsetY = 0;
	if (sqrt((aop6[0] * aop6[0] + aop6[1] * aop6[1] + aop6[3] * aop6[3] + aop6[4] * aop6[4]) / 2) < 0.8) {
		offsetX = -0.5; offsetY = -0.5;
	}
	switch (mode)
	{
	case BILINEAR:
		for (r = 0; r < nBlockRows; r++){
			if (stBlockRow + r>nReRows - 1) break;
			//		x2 = x1;	y2 = y1;
			x1 = r*aop6[1]; y1 = r*aop6[4];
			for (c = 0; c < nBlockCols; c++/*, x2 += dxx, y2 += dyx*/){
				if (stBlockCol + c>nReCols - 1) break;
				x2 = c*aop6[0] + x1 + aop6[2] + offsetX;
				y2 = c*aop6[3] + y1 + aop6[5] + offsetY;
				if (x2 < 0 || x2 >= nCols - 1 || y2 < 0 || y2 >= nRows - 1){
					pBuf[c] = 0;
					continue;
				}
				nx = int(x2);	ny = int(y2);
				float dx = x2 - nx;
				float dix = 1 - dx;
				float dy = y2 - ny;
				float diy = 1 - dy;
				const BYTE* pData = data + nCols*ny + nx;
				float xval1 = pData[0] * dix + pData[1] * dx;	pData += nCols;
				float xval2 = pData[0] * dix + pData[1] * dx;
				pBuf[c] = BYTE(xval1*diy + xval2*dy);
			}
			pBuf += nReCols;
			//		x1 += dxy;	y1 += dyy;
		}
		break;
	case CUBIC:
		for (r = 0; r < nBlockRows; r++){
			if (stBlockRow + r>nReRows - 1) break;
			//		x2 = x1;	y2 = y1;
			x1 = r*aop6[1]; y1 = r*aop6[4];
			for (c = 0; c < nBlockCols; c++/*, x2 += dxx, y2 += dyx*/){
				if (stBlockCol + c>nReCols - 1) break;
				x2 = c*aop6[0] + x1 + aop6[2] + offsetX;
				y2 = c*aop6[3] + y1 + aop6[5] + offsetY;
				if (x2 < 1 || x2 >= nCols - 2 || y2 < 1 || y2 >= nRows - 2){
					pBuf[c] = 0;
					continue;
				}
				nx = int(x2);	ny = int(y2);
				float dx = x2 - nx;	float dx_2 = dx*dx; float dx_3 = dx_2*dx;				
				float dy = y2 - ny;	float dy_2 = dy*dy; float dy_3 = dy_2*dy;
				float wx1 = -dx + 2 * dx_2 - dx_3;
				float wx2 = 1 - 2 * dx_2 + dx_3;
				float wx3 = dx + dx_2 - dx_3;
				float wx4 = -dx_2 + dx_3;
				float wy1 = -dy + 2 * dy_2 - dy_3;
				float wy2 = 1 - 2 * dy_2 + dy_3;
				float wy3 = dy + dy_2 - dy_3;
				float wy4 = -dy_2 + dy_3;
				const BYTE* pData = data + nCols*(ny-1) + nx - 1;
				float xval1 = pData[0] * wx1 + pData[1] * wx2 + pData[2] * wx3 + pData[3] * wx4;	pData += nCols;
				float xval2 = pData[0] * wx1 + pData[1] * wx2 + pData[2] * wx3 + pData[3] * wx4;	pData += nCols;
				float xval3 = pData[0] * wx1 + pData[1] * wx2 + pData[2] * wx3 + pData[3] * wx4;	pData += nCols;
				float xval4 = pData[0] * wx1 + pData[1] * wx2 + pData[2] * wx3 + pData[3] * wx4;	
				wx1 = xval1*wy1 + xval2*wy2 + xval3*wy3 + xval4*wy4;
				if (wx1 < 0) pBuf[c] = 0; else if (wx1>255) pBuf[c] = 255;
				else pBuf[c] = BYTE(wx1);				
			}
			pBuf += nReCols;
			//		x1 += dxy;	y1 += dyy;
		}
		break;
	}
	
	return true;
}

bool Pretreatment::WriteBlock2Buf(BYTE* pBuf, int nCols, int nRows, int stBlockCol, int stBlockRow,const BYTE* pBlock, int nBlockCols, int nBlockRows)
{
	if (stBlockCol + nBlockCols > nCols) nBlockCols = nCols - stBlockCol;
	if (stBlockRow + nBlockRows > nRows) nBlockRows = nRows - stBlockRow;

	pBuf = pBuf + stBlockRow*nCols + stBlockCol;
	for (int r = 0; r < nBlockRows; r++){
		memcpy(pBuf, pBlock, sizeof(BYTE)*nBlockCols);
		pBuf	+= nCols;
		pBlock += nBlockCols;
	}
	return true;
}

bool Pretreatment::ReadBlock4Buf(const BYTE* pBuf, int nCols, int nRows, int stBlockCol, int stBlockRow, BYTE* pBlock, int nBlockCols, int nBlockRows)
{
	if (stBlockCol<0) stBlockCol = 0;
	if (stBlockRow<0) stBlockRow = 0;
	if (stBlockCol + nBlockCols > nCols) nBlockCols = nCols - stBlockCol;
	if (stBlockRow + nBlockRows > nRows) nBlockRows = nRows - stBlockRow;
	pBuf = pBuf + stBlockRow*nCols + stBlockCol;
	for (int r = 0; r < nBlockRows; r++){
		memcpy(pBlock, pBuf, sizeof(BYTE)*nBlockCols);
		pBuf += nCols;
		pBlock += nBlockCols;
	}
	return true;
}

inline int get_split_num(double range, double blockSz,double minBlockRatio = 0.1 )
{
	if (range<0) return 0;
	double	fSplit = range / blockSz;
	int		nSplit = int(fSplit);
	if (fSplit - nSplit > minBlockRatio) nSplit++;
	return nSplit;
}

double* Pretreatment::split_rectangle(int& nColSplitNum, int& nRowSplitNum, double rect[4], double max_nCols, double max_nRows, double overlay_ratio /* = 0 */)
{
	if (overlay_ratio >= 1 || overlay_ratio < 0) return NULL;

	double xmax, xmin, ymax, ymin;
	xmin = rect[0];	ymin = rect[1];
	xmax = rect[2];	ymax = rect[3];

	double dColRange = xmax-xmin;
	double dRowRange = ymax-ymin;

	double nCols_cover		= max_nCols*overlay_ratio;	double nRows_cover = max_nRows*overlay_ratio;
	double nCols_nocover = max_nCols - nCols_cover;	double nRows_nocover = max_nRows - nRows_cover;

	nColSplitNum = dColRange <= max_nCols ? 1 : get_split_num(dColRange - nCols_cover, nCols_nocover);
	nRowSplitNum = dRowRange <= max_nRows ? 1 : get_split_num(dRowRange - nRows_cover, nRows_nocover);

	double* recBuf = new double[4 * nColSplitNum*nRowSplitNum];

	int r, c;
	double stCol, stRow;

	double* pBuf = recBuf;
	for (r = 0, stRow = ymin; r < nRowSplitNum; r++, stRow += nRows_nocover){
		stCol = xmin;
		for (c = 0; c < nColSplitNum; c++){
			*pBuf = stCol;	pBuf++;
			*pBuf = stRow;	pBuf++;
			*pBuf = stCol + max_nCols<=xmax ? max_nCols : xmax - stCol;	pBuf++;
			*pBuf = stRow + max_nRows<=ymax ? max_nRows : ymax - stRow;	pBuf++;
			stCol += nCols_nocover;
		}
	}

	return recBuf;
}

double* Pretreatment::split_rectangle(double rect[4], int nColSplitNum, int nRowSplitNum, double overlay_ratio /* = 0 */)
{
	if (overlay_ratio >= 1 || overlay_ratio < 0) return NULL;

	double xmax, xmin, ymax, ymin;
	xmin = rect[0];	ymin = rect[1];
	xmax = rect[2];	ymax = rect[3];

	double range_x = (xmax - xmin) / (nColSplitNum + overlay_ratio);
	double range_y = (ymax - ymin) / (nRowSplitNum + overlay_ratio);
	double x_cover = range_x*overlay_ratio;	double y_cover = range_y*overlay_ratio;
	double x_nocover = range_x - x_cover;	double y_nocover = range_y - y_cover;

	double* recBuf = new double[4 * nColSplitNum*nRowSplitNum];

	int r, c;
	double stx, sty;

	double* pBuf = recBuf;
	for (r = 0, sty = ymin; r < nRowSplitNum; r++, sty += y_nocover){
		stx = xmin;
		for (c = 0; c < nColSplitNum; c++){
			*pBuf = stx;	pBuf++;
			*pBuf = sty;	pBuf++;
			*pBuf = range_x;	pBuf++;
			*pBuf = range_y;	pBuf++;
			stx += x_nocover;
		}
	}

	return recBuf;
}

int* Pretreatment::split_image(int& nColSplitNum, int& nRowSplitNum, int stCol, int stRow, int nCols, int nRows, int nBits, int maxBlockSz, SPLITMODE mode /* = COL_MAX_FIRST */, double overlay_ratio /* = 0 */)
{
	int nBytes = nBits / 8;		int nBlockSz = maxBlockSz / nBytes;
	int minBlockEdgeLen = 100;

	float dataSz = nRows*1.0f*nCols;
	double numratio =  dataSz / nBlockSz;
	int blockNum;

	if (numratio <= 1 ) { 
		if (mode != RECT_FIRST){
			lp1:
			int* block = new int[4];
			block[0] = stCol;	block[1] = stRow;
			block[2] = nCols;	block[3] = nRows;
			nColSplitNum = nRowSplitNum = 1;
			return block;
		}
		else{
			float whratio = (float)nCols / nRows;
			if (whratio > 1 ){
				nColSplitNum = int(whratio + 0.5);
				nRowSplitNum = 1;
				goto lp2;
			}else if (whratio<1){
				nColSplitNum = 1;
				nRowSplitNum = int(1/whratio+0.5);
				goto lp2;
			}
			else {
				nColSplitNum = nRowSplitNum = 1;
				goto lp1;
			}
		}
	}

	blockNum = int(numratio);	
	if (numratio-blockNum > 0.00001) blockNum++;
//	if (numratio - blockNum>0.1) blockNum++;

//	int nBlockCols, nBlockRows;
	if (mode == RECT_FIRST){
		int len = int(sqrt((float)nBlockSz));
		nColSplitNum = nCols / len;	if (nColSplitNum < 1) nColSplitNum = 1;
		nRowSplitNum = nRows / len;	if (nRowSplitNum < 1) nRowSplitNum = 1;
		while (nColSplitNum*nRowSplitNum*1.0f*nBlockSz < dataSz)
		{
			if (len>nCols) nRowSplitNum++;
			else if (len>nRows) nColSplitNum++;
			else{
				float r = nCols*1.0f*nRowSplitNum / (nColSplitNum * nRows*1.0f ) ;
				if (r>1.2f) nColSplitNum++;
				else if(r<0.833f) nRowSplitNum++;
				else {
					nColSplitNum++;
					nRowSplitNum++;
				}
			}
		}
	}
	else if (mode==WIDTH_HEIGH_RATIO_FIRST){
		numratio = sqrt(numratio);
		int tmp = int(numratio);	if (numratio - tmp > 0.00001) tmp++;
		nColSplitNum = nRowSplitNum = tmp;
	}
	else
	{
		int nCol0, nRow0;
		int *pNumC, *pNumR;
		if (mode==COL_MAX_FIRST){
			nCol0 = nCols;
			nRow0 = nRows;
			pNumC = &nColSplitNum;
			pNumR = &nRowSplitNum;
		}
		else
		{
			nCol0 = nRows;
			nRow0 = nCols;
			pNumC = &nRowSplitNum;
			pNumR = &nColSplitNum;
		}
		if (minBlockEdgeLen*nCol0 < nBlockSz) {
			*pNumC = 1;
			*pNumR = blockNum;
		}
		else
		{
			double tmp = (nRow0 - overlay_ratio*minBlockEdgeLen) / (minBlockEdgeLen - overlay_ratio*minBlockEdgeLen);
			*pNumR = int(tmp);	if (tmp - *pNumR>0.00001) *pNumR++;
			tmp = numratio / *pNumR;
			*pNumC = int(tmp);	if (tmp - *pNumC > 0.00001) *pNumC++;
		}
	}
// 	nBlockRows = nBlockSz / nCols;	if (nBlockRows < minBlockRows) nBlockRows = minBlockRows;
// 	if(nBlockRows<nRows) nBlockRows = nRows / (nRows / nBlockRows);
// 	nBlockCols = nBlockSz / nBlockRows;
// 
 	
// 	double* rect = split_rectangle(nColSplitNum, nRowSplitNum, rc, nBlockCols, nBlockRows, overlay_ratio);
lp2:
	double rc[4] = { stCol, stRow, stCol + nCols, stRow + nRows };
	double* rect = split_rectangle(rc,nColSplitNum, nRowSplitNum, overlay_ratio);
	if (!rect) return NULL;

	int nRectNum = nColSplitNum*nRowSplitNum;
	int* block = new int[4 * nRectNum];

	double* pRect = rect;	int* pBlock = block;
	for (int i = 0; i < 4*nRectNum; i++){
		*pBlock = int(*pRect);
		pRect++;	pBlock++;
	}
	delete rect;
	return block;
}

int* Pretreatment::split_image(int stCol,int stRow, int nCols, int nRows, int nColSplitNum, int nRowSplitNum, double overlay_ratio /* = 0 */)
{
	double rc[4] = { stCol, stRow, stCol+nCols, stRow+nRows };
	double* rect = split_rectangle(rc,nColSplitNum,nRowSplitNum,overlay_ratio);
	if (!rect) return NULL;

	int nRectNum = nColSplitNum*nRowSplitNum;
	int* block = new int[4 * nRectNum];
	
	double* pRect = rect;	int* pBlock = block;
	for (int i = 0; i < 4*nRectNum; i++){
		*pBlock = int(*pRect);
		pRect++;	pBlock++;
	}
	delete rect;
	return block;
}