#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
from os.path import isfile, join

import cv2


def video(pathIn, pathOut, fps):
    '''
    Generates output video with given paths
    pathIn  - path of video frames
    pathOut - path of video and name
    fps     - rate of frames per second
    '''
    frame_array = []
    files = [f for f in os.listdir(pathIn) if isfile(join(pathIn, f))]
    #for sorting the file names properly
    files.sort()
    
    for i in range(len(files)):
        filename = pathIn + '/' + files[i]
        
        #reading each files
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        
        #inserting the frames into an image array
        frame_array.append(img)
        
#    fourcc = VideoWriter::fourcc('p', 'n', 'g', ' ')
#    fourcc = cv2.VideoWriter_fourcc('F','F','V','1')
#    fourcc = cv2.VideoWriter_fourcc('M','J','P','G')
    fourcc = 0
    out = cv2.VideoWriter(pathOut, fourcc, fps, size)
    
    for i in range(len(frame_array)):
        # writing to a image array
        if i == len(frame_array)-1:
            out.write(frame_array[i])
            out.write(frame_array[i])
            out.write(frame_array[i])
            out.write(frame_array[i])
        else:
            out.write(frame_array[i])
        
    out.release()


def video_output_nonstat(xds_out, p_figs, fps):
    '''
    Generates paths for input frames and output video
    xds_out - output dataset
    p_figs  - path of animation video frames
    fps     - rate of frames per second
    '''
    for case_ix in xds_out.case.values[:]:

        # select output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        name = 'output_video_{0}.avi'.format(case_id)
        pathIn= op.join(p_figs, case_id)
        pathOut = op.join(p_figs, name)

        video(pathIn, pathOut, fps)
    

def video_vortex(xds_out, p_figs, fps):
    '''
    Generates paths for input frames and output video
    p_figs  - path of animation video frames
    fps     - rate of frames per second
    '''
#    for case_ix in xds_out.case.values[:]:

        # select output case subfolder
#        case_id = '{0:04d}'.format(case_ix)
#        name = 'vortex_video_{0}.avi'.format(case_id)
#        pathIn= op.join(p_figs, case_id, 'vortex_figs', '0000')
    name = 'vortex_video_0.avi'
    pathIn= op.join(p_figs, 'output_figs_vortex')
    pathOut = op.join(p_figs, name)

    video(pathIn, pathOut, fps)
#    pathIn= op.join(p_figs, '0000')
#    pathOut = op.join(p_figs, 'vortex_video.avi')


def video_output_panel(xds_out, p_figs, fps):
    '''
    Generates paths for input frames and output video
    xds_out - output dataset
    p_figs  - path of animation video frames
    fps     - rate of frames per second
    '''
    for case_ix in xds_out[0].case.values[:]:

        # select output case subfolder
        case_id = '{0:04d}'.format(case_ix)
        name = 'panel_video_{0}.avi'.format(case_id)
        pathIn= op.join(p_figs, case_id)
        pathOut = op.join(p_figs, name)

        video(pathIn, pathOut, fps)


