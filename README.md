# Self-Pinning

This folder contains the Matlab scripts and resulting data used to create the figures in the publication

*Self-Pinning Transition of a Tonks-Girardeau Gas in a Bose-Einstein Condensate*  
Tim Keller, Thomás Fogarty, and Thomas Busch  
Physical Review Letters **128**, 053401 (2022).

Simulations were performed with *MATLAB R2019a* and can be executed e.g. via

	cd(fullfile(pwd,'/matlab/'))
	run('selfpinning_Fig1.m')
	run('selfpinning_Fig2.m')
	run('selfpinning_Fig3.m')
	run('selfpinning_FigS1.m')
	run('selfpinning_FigS2.m')
	run('selfpinning_FigS3.m')

External dependencies are the additional files *'fftdef.m'* and *'v2struct.m'* which are included in the *matlab/* folder. 
The colormap used in Fig. 2 is [*'Smooth Cool Warm'*](https://www.kennethmoreland.com/color-advice/) by Kenneth Moreland.

The figures are created via tikz in the .tex file for the publication itself. 
Tikz settings need to be included in the preamble via

	\input{figures/tikz_settings.tex}

The figures are then created at the desired locations via

	\begin{figure}
	\centering
	\input{figures/selfpinning_Fig1.tikz}
	\caption{Caption}
	\label{fig:Fig1}
	\end{figure}
	
For questions please contact tim.keller@oist.jp 	
