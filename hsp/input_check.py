from flask import Flask, request, jsonify
import os
import random
from osgeo import gdal

ALLOWED_EXTENSIONS = {'tif'}


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def get_image_files(folder_path):
    if not os.path.isdir(folder_path):
        return None

    image_files = []
    for filename in os.listdir(folder_path):
        if allowed_file(filename):
            image_files.append(filename)
    image_files.sort()
    return image_files


def analyze_images(folder_path):
    image_files = get_image_files(folder_path)
    if not image_files:
        return None

    total_count = len(image_files)

    sorted_files = sorted(image_files)
    damaged_files = []
    damaged_reasons = []
    for filename in sorted_files:
        file_path = os.path.join(folder_path, filename)
        ds = gdal.Open(file_path)
        try:
            x = ds.RasterXSize
        except:
            damaged_files.append(filename)
            damaged_reasons.append("file damaged")
            continue
        basename = os.path.basename(file_path).split('.')[0]
        rpb_name = basename + '.rpb'
        rpb_file = os.path.join(folder_path, rpb_name)
        if not os.path.exists(rpb_file):
            damaged_files.append(filename)
            damaged_reasons.append("lack rpb")


    return {
        'total_count': total_count,
        'images_list': sorted_files,
        'damaged_files': damaged_files,
        'damaged_reasons':damaged_reasons,
        'damaged_count': len(damaged_files),
        'valid_count': total_count - len(damaged_files)
    }



