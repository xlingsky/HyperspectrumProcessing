import cv2 as cv
import numpy as np
from tqdm import tqdm

def computeRowShifts(image, window_size=31, inner_window=15, method = cv.TM_CCOEFF_NORMED, threshold_response = 0.8 , threshold_offset_std = 3):
    H, W = image.shape

    offsets = np.zeros(H)
    offsets_std = np.zeros(H)
    inner_st = (window_size-inner_window) // 2
    for i in tqdm(range(1, H), desc="Computing offsets"):
        last_row = image[i - 1, :]
        cur_row = image[i, :]
        row_offsets = np.zeros(W//window_size)
        row_confidence = np.zeros(W//window_size)
        cnt = 0
        for j in range(window_size, W, window_size):
            st = j - window_size
            source = last_row[st:j]
            template = cur_row[st+inner_st:st+inner_st+inner_window]
            res = cv.matchTemplate( source, template, method)
            minv, maxv, min_loc, max_loc = cv.minMaxLoc(res)
            if method in [cv.TM_SQDIFF, cv.TM_SQDIFF_NORMED]:
                row_offsets[cnt] = min_loc[1]
                row_confidence[cnt] = minv
            else:
                row_offsets[cnt] = max_loc[1]
                row_confidence[cnt] = maxv

            row_offsets[cnt] -= inner_st
            cnt += 1

        row_offsets_valid = row_offsets[row_confidence > threshold_response]
        if len(row_offsets_valid) == 0: continue
        offsets_std[i] = np.std(row_offsets_valid)
        offsets[i] = np.mean(row_offsets[:cnt])
    offsets[offsets_std > threshold_offset_std] = 0
    offsets = -np.round(offsets).astype(np.int32)
    return offsets

def computeRowShifts_pc(image, threshold = 0.8):
    H, W = image.shape
    row = 32
    col = W//row
    l = row * col

    hann = cv.createHanningWindow((col, row), cv.CV_32F)

    shifts = np.zeros(H)
    responses = np.zeros(H)
    for i in tqdm(range(1, H), desc="Computing offsets"):
        last_row = image[i - 1, :l].reshape((row, col)).astype(np.float32)
        cur_row = image[i, :l].reshape((row, col)).astype(np.float32)

        shift, res = cv.phaseCorrelate(cur_row, last_row, hann)
        shifts[i] = shift[0]
        responses[i] = res

        # fft_sig1 = np.fft.fft(last_row)
        # fft_sig2 = np.fft.fft(cur_row)
        # fft_sig2_conj = np.conj(fft_sig2)

        # R = (fft_sig1 * fft_sig2_conj) / abs(fft_sig1 * fft_sig2_conj)
        # r = np.fft.ifft(R)

        # shift = np.argmax(r)
        # response = np.max(r)
        # if response > threshold:
        #     offsets[i] = shift

    shifts[responses < threshold] = 0
    shifts = np.round(shifts).astype(np.int32)
    return shifts

def rectify(image, shifts, width = None):
    H, W = image.shape
    if width is None: 
        width = W - min(shifts)+max(shifts)
        shifts = shifts - min(shifts)

    rectified = np.zeros((H, width), dtype=image.dtype)
    for i in range(H):
        rectified[i, shifts[i]:min(width, W+shifts[i])] = image[i, :min(width-shifts[i], W)]

    return rectified

if __name__ == "__main__":
    import sys, os
    imgpath = sys.argv[1]

    img = cv.imread(imgpath, cv.IMREAD_GRAYSCALE)
    
    shifts_relative = computeRowShifts(img)
    shifts_absolute = np.zeros_like(shifts_relative)
    for i in range(1, len(shifts_relative)):
        shifts_absolute[i] = shifts_absolute[i-1] + shifts_relative[i]

    rectified = rectify(img, shifts_absolute)

    cv.imwrite(os.path.splitext(imgpath)[0] + "_rectified.png", rectified)
