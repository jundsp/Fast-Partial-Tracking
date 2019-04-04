# Fast Partial Tracking

Matlab code that fully implements the partial tracking method published in the following paper:

J. Neri and P. Depalle, "<a href="http://dafx2018.web.ua.pt/papers/DAFx2018_paper_26.pdf" target="_blank">**Fast Partial Tracking of Audio with Real-Time Capability through Linear Programming**</a>," In *Proceedings of the 21st International Conference on Digital Audio Effects (DAFx-18)*, Aveiro, Portugal, pp. 326-333, Sep. 2018.

The script `jun_demo.m` runs a demonstration of the partial tracking method, followed by a re-synthesis of the signal from the detected partials. Instantaneous parameters of the partials are written to a binary file after each analysis frame, which can be read into Matlab again with the corresponding read function.

<img src="https://www.music.mcgill.ca/~julian/dafx18/pictures/partials_kara.png" alt="" width="500">
