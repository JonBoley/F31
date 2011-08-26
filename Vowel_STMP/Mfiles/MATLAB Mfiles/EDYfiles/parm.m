Content of vector parm. Values given as 0 are not used.

  1 - 0
  2 - 0 if 3 lines of stim. spec. text (L. M. R.)
      ‚0 if short-form, 1 line only (S.)
  3 - 0
  4 - Name of program that took the data: 'ST', 'TC', 'AV', etc.
  5 - 0
  6:10 - date as '10.28.97'
 11:17 - track, unit, picture no. ' 3.05 ,P47    '
 18 - 0
 19:38 - L.or S. stimulus specification text

STIMULUS PARAMETERS
 39 - left stimulus: 0-tone, 3-noise, 4-click, 6-RMF,
      7-Tucker-Davis, 8-ROM farm, 9-noisebands, 10-2bands,
	  11-no stimulus
 40 - left stimulus: PF tone freq., -1 freq. range (see 69:73),
      ROM or RMF sync freq (Hz), =nn for 'XRnn', otherwise 0
 41:42 - same for right stimulus

 43:44 - left attenuation: (low atten, high atten) or (atten, -1)
 45:46 - same for right attenuation
 
 47:50 - left timing: (duration ms, PF low delay (x10 ms),
         PF high delay or 0 if no range, gate (-1 NG, 0 UG, 1 gated)
 51:54 - same for right timing
 55 - 0 or attenuation step size in TC

FILE DESCRIPTION
 56:60 - file name for this picture 'A12345.DAT'
 61 - number of blocks (512 bytes) in the data
 62 - picture number
 63 - 1
 64 - data format:   0 end of experiment       1 dot display (ST)
                     2 TD format (not used)    3 rate-level (not used)
					 4 tuning curve           15 at exit from SE
					19 analog averager AVELN) 30 acoustic calib (SRCAL)
					31 JJR acoustic data     „64 Ryugo lab files
 65:68 - 0

MORE STIMULUS PARAMETERS
 69:73 - left ch. frequency range: PF low freq, PF high freq, PF center
         freq, step size (Hz if linear, (mult-1)*32768 if log), 1 if log
		 or 0 if linear
 74:78 - same for right channel
 
 79:83 - left ch. ROM or RMF parameters; start block, number of blocks,
         10(ROM) or 1000(RMF)*syncfreq (Hz), rate divisor, 0 or RMF channel
 84:88 - same for right channel
 
 89 - spectral level, in dB re tone level, of noise stim in left ch.
 90 - same for right ch.
 
MODE-LINE AND GROUP PARAMETERS
 91 - number of repetitions of each stimulus
 92 - ear: 1-RI, 2-LE, 3-LE+RI, 4-ROT, 5-RIM, 6-LEM, 7-BO
 93 - presentation order 0-SEQ, 1-RND, 2-CEN
 94 - stimulus repetition period, ms
 95 - sync source: 0-no sync, 1-SYNR, 2-SYNL
 96 - 0  (1 if EP recorded during spike data)
 97 - type of stimulus range:   0-no range     1-frequency     2-atten.
                                3-delay        4-ear (ROT)     5-Tucker/
								                           Davis file rotation
	   (TDT noise band freq sweeps coded as 1)
 98 - section count:  1 (RND, CEN)  or no. of reps/stimulus (SEQ)

PARAMETERS MEANINGFUL IN STIM ONLY
 99:108 - stimulus control words
 
109 - data resolution for spike timing, microseconds/clock tick

TC PARAMETERS
110 - number of data points in the TC
111 - number of trials meeting criterion for threshold
112 - ratio of search frequency step size to TC frequency step size
113 - index of BF in the TC data
114:115 - threshold (real*4) dB
116:117 - BF (real*4)
118:119 - spontaneous rate (real*4)
120:121 - total number of stimuli (real*4)
122 - 10*spont rate
123 - criterion, driven - spont spikes for a threshold

124:128 - 0

MODE AND RIGHT CHANNEL STIMULUS SPECIFICATION TEXT (only if parm(2)=0)
129:148 - mode line
149:168 - right line
169 - 0

TEXT LINES FOR AV AND FOR PC MESSAGE FOR Tucker/Davis STIMULI
170:182 - AV text line (not used in current version AVELN)
170:188 - PC communication text line; created in SE, sent by ST to
          PC programs (commains or filrot)

189:240 - 0  (used for parameters in Ryugo lab system)

MISC PARAMETERS
241:242 - PC's Tucker-Davis random number generator kernal
243:244 - 0 (used for params in Ryugo system)
245:249 - electronic switch and digital oscillator control words
250:251 - SE's random number generator kernal

252:255 - time of day at start of picture  'hh:mm:ss'

256 - 0
