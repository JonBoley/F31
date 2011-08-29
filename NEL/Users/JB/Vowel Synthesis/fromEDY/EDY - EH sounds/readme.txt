     This disk contains stimulus data files for 6 synthetic speech sounds.
The data files are unformatted, direct access files.  The data are
written as 16 bit integers in twos-complement format, EXCEPT THAT BIT 0
IS NOT A DATA BIT.  Bit 0 is a sync mark which occurs once in each repetition
period of the stimulus.

     The stimuli and their characteristics are listed below:

     File    Stimulus       Fund.      Sync        Sampling      CCD
                            freq.      freq.         freq.       word
   --------------------------------------------------------------------
    EH1.DAT  vowel EH       112 Hz     112 Hz       57604 Hz    036007
             as in /met/

    EH2.DAT  vowel EH       128 Hz     128 Hz       32468 Hz   0114005
             as in /met/

    AH.DAT   vowel AH       128 Hz     128 Hz       32468 Hz   0114005
             as in /father/

    IH.DAT   vowel IH       128 Hz     128 Hz       32468 Hz   0114005
             as in /bit/

    BA.DAT   CV /ba/        not        100 Hz        9968 Hz   0161013
                            periodic

    DA.DAT   CV /da/        not        100 Hz        9968 Hz   0161013
                            periodic

The lengths of the files varies with the length of the stimuli.  Only 
one or two periods of the periodic stimuli are contained in their files.
These should be concatenated to make a continuous periodic stimulus.  
The /ba/ and /da/ stimulus files contain the whole stimulus (approx. 
100 mS long plus some zero padding).


     The files VOW*.DAT are the stimuli used by Conley and Keilson.  They
are /EH/-like vowels with the following parameters:

     F0 = 100.16 Hz		F1 = 500.8 Hz (BW 60 Hz)
     F2 = varies (90 Hz)	F3 = 2500 Hz (200 Hz)
     F4 = 3300 Hz (250 Hz)	F5 = 3750 Hz (200 Hz)
(done with ACSL version of Klatt, see EDY's vowel data binder for details)

     File    Stimulus       Fund.      Sync        Sampling       F2
                            freq.      freq.         freq.       freq.
   --------------------------------------------------------------------
   VOW14		   100.16     no syncs     51.282        1.4 kHz
   VOW15		     "           "           "           1.5 kHz
   VOW16		     "           "           "           1.6 kHz
   VOW17      /EH/           "           "           "           1.7 kHz
   VOW18                     "           "           "           1.8 kHz
   VOW19		     "           "           "           1.9 kHz
   VOW20		     "           "           "           2.0 kHz

In 2/94 additional vowles were synthesized with same parameters, except
F2, these have the following F2s:
   1450, 1550, 1625, 1650, 1675, 1700, 1725, 1750, 1775, 1850, 1950
   2050, 2100, 2150, and 2200 Hz.
These were written to the Behavioral Chamber ROMs on 3/15/94 as below.
On 6/14/94 three additional vowels in this series were synthesized, with
F2=1701, 1705, and 1710 Hz.  They were programmmed on 6/15/94.
They are in the files VOnnnn.DAT where nnnn = F2.  In the Suns, they are
call BOBHEIHZ vowels.

The files BEHRHI.DAT and BEHRLO.DAT contain the upper and lower bytes
of behavioral chamber ROMs as of 3/15/94.  These contain the following:

    Vowel     start      start		where start block means block
              block      address	number in decimal given to SELECT
  ---------------------------------	as in ROM xx/2/195, block number
    /EH/	 0	     0		is xx.
    /ae/	 2	 01000
    /AH/	 4	 02000		and start address is start address
    /AW/	 6	 03000		of 01000 words in vowel in ROM
    /UH/	 8	 04000		programmer program.  Equals starting
  F2=2.0 kHz	10	 05000		address on ROM in octal.
    =1.9 kHz	12	 06000
    =1.8 kHz	14	 07000
    =1.7 kHz	16	010000
    =1.65 kHz   18	011000   (new)
    =1.6 kHz	20	012000
    =1.55 kHz 	22	013000   (new)
    =1.5 kHz    24	014000
    =1.45 kHz	26	015000   (new)
    =1.4 kHz	28	016000
    =1.625 kHz	30	017000   (ALL NEW from here on)
    =1.675 kHz  32      020000
    =1.7 kHz    34      021000   (differs from 1.7 kHz above by time delay)
    =1.725 kHz  36      022000
    =1.75 kHz   38      023000
    =1.85 kHz   40      024000
    =1.95 kHz   42      025000
    =2.0  kHz   44      026000   (differs from 2.0 kHz above by time delay)
    =2.05 kHz   46      027000
    =2.1 kHz    48      030000
    =2.15 kHz   50      031000
    =2.2 kHz    52      032000
    =1.701 kHz  54      033000
    =1.705 kHz  56      034000
    =1.71 kHz   58      035000

The parameters of the vowels in blocks 10-58 are given above (VOWxx.DAT).
The parameters of the first 5 vowels are unknown.
