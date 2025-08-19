#include "bigfileio.h"


#if defined(WIN32)||defined(_WIN32)
HANDLE CreateFileE(LPCSTR lpstrPathName, UINT dwDesiredAccess){
	if (dwDesiredAccess == GENERIC_READ)
		return ::CreateFile(lpstrPathName, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING, FILE_FLAG_RANDOM_ACCESS, NULL);
	else return ::CreateFile(lpstrPathName, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, 0, NULL);
}
#else

void CloseHandle(int hFile){ close(hFile); }

int WriteFile(int hFile, const void *pBuf, DWORD bufSize, DWORD *opSz, int *lpOverlapped){
	*opSz = write(hFile, pBuf, bufSize);
  return *opSz>0?1:0;
}

int ReadFile(int hFile, void *pBuf, DWORD bufSize, DWORD *opSz, unsigned int *lpOverlapped){
	*opSz = read(hFile, pBuf, bufSize);
  return *opSz>0?1:0;
}

unsigned int SetFilePointer(int hFile, long lDistanceToMove, long *lpDistanceToMoveHigh, int dwMoveMethod){
	OFFSET ret, dis = 0; if (lpDistanceToMoveHigh) dis = *lpDistanceToMoveHigh;
	dis = (dis << 32) | lDistanceToMove;
	ret = lseek64(hFile, dis, dwMoveMethod);
	return (unsigned int)ret;
}

HANDLE CreateFileE(LPCSTR lpstrPathName, UINT dwDesiredAccess){
	if (dwDesiredAccess == GENERIC_READ) return open64(lpstrPathName, O_RDWR);
	return open64(lpstrPathName, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
}

#endif
