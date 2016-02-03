Blocks
======

Present-day deformation at plate boundary zones, recorded as GPS velocities, reflect the combination of plate motion and elastic strain accumulation. This suite of codes implements the block modeling methodology described by Meade and Loveless (2009), Block modeling with connected fault-network geometries and a linear elastic coupling estimator in spherical coordinates, Bulletin of the Seismological Society of America.

Blocks is designed for use with Matlab R2014b and later. 

Documentation (evolving): https://docs.google.com/document/d/1AJheJrVqPX4yj2hbgysC-H2RkdkfkxvNJQIza1b3u34/edit?usp=sharing

To get started, run the following commands on the Matlab command prompt:

    cd('~/MATLAB/Blocks') % Edit path to Blocks directory
    cd BlocksUtilities
    blockspath % This function adds the Blocks subdirectories to your Matlab path

You can create a new template model directory structure using:

    blocksdirs('~/MATLAB/Blocks/California') 
    % Edit path to your project name; a new directory will be created if it doesn't exist

Then, edit the Blocks geometry files (.segment and .block) using SegmentManager:

    cd ~/MATLAB/Blocks/California/command
    SegmentManager
    % Within SegmentManager, click "Load" under "Command file" and load 'model.command'. 
    % Use SegmentManager tools to add and modify segment and block properties, saving 
    % the geometry files to the ../segment and ../block directories
To run the analysis,

    cd ../result
    Blocks('../command/model.command')
    % The results will be saved in a newly generated directory in the result directory
To view the results,

    ResultManager
    % Load a result directory. If you have more than one set of results, you can compare
    % them by loading both a "Result directory" and "Compare directory"

SegmentManager interface:
-------------------------
![alt tag](https://cloud.githubusercontent.com/assets/4225359/9386297/d46874ca-4728-11e5-9deb-48899bd91770.png)
