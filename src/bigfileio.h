#ifndef BIGFILEIO_H
#define BIGFILEIO_H

#if defined(WIN32)||defined(_WIN32)

#include "windows.h"

#ifndef _CreateFileE
#define _CreateFileE
HANDLE CreateFileE(LPCSTR lpstrPathName, UINT dwDesiredAccess);
#endif
#else
#define INVALID_HANDLE_VALUE			 (HANDLE)-1

typedef 	unsigned int	HANDLE		;
typedef 	int				BOOL		;
typedef 	const char*     LPCSTR		;
typedef		char*			LPSTR		;
typedef 	unsigned char	BYTE		;
typedef 	unsigned int	RGBQUAD		;
typedef 	unsigned short	WORD		;
typedef 	unsigned int	DWORD		;
typedef 	unsigned int	UINT		;

#define __USE_LARGEFILE64
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>

#ifdef __APPLE__
#define lseek64 lseek
#define open64   open
#endif

typedef off_t OFFSET;

typedef union _LARGE_INTEGER {
	struct { DWORD LowPart; unsigned int HighPart; }u;
	long QuadPart;
}LARGE_INTEGER;

#define GENERIC_READ                     (0x80000000L)
#define GENERIC_WRITE                    (0x40000000L)
#define FILE_BEGIN						 SEEK_SET
#define FILE_CURRENT					 SEEK_CUR
#define FILE_END						 SEEK_END
#define HFILE_ERROR						 (-1)

void CloseHandle(int hFile);
int WriteFile(int hFile, const void *pBuf, DWORD bufSize, DWORD* opSz, int *lpOverlapped);
int ReadFile(int hFile, void *pBuf, DWORD bufSize, DWORD* opSz, unsigned int *lpOverlapped);
unsigned int SetFilePointer(int hFile, long lDistanceToMove, long *lpDistanceToMoveHigh, int dwMoveMethod);
HANDLE CreateFileE(LPCSTR lpstrPathName, UINT dwDesiredAccess);

#endif // end for WIN32

#endif
