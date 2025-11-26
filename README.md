
<h1 align="center"> On the relation of skin stretch and finger joint angle evolution in human hand grasping tasks</h1>

 <p align="center">

  <img width="400" src="https://github.com/EleFontana/SkinStretchEval/blob/main/media/cp_logo.png" alt="cp_logo">
  <img width="300" src="https://github.com/EleFontana/SkinStretchEval/blob/main/media/unipi_logo.png" alt="unipi_logo_nero">
  <img width="300" src="https://github.com/EleFontana/SkinStretchEval/blob/main/media/Logo_DII.jpg" alt="dii_logo"> 
  <img width="300" src="https://github.com/EleFontana/SkinStretchEval/blob/main/media/iit_logo.jpg" alt="iit_logo"> 
</p>

# Abstract
There is strong evidence that skin stretch at joint level contributes to proprioception, i.e., the inner perceptual representation of body position and motion. Considering the dorsum of the hand and fingers, state-of-the art studies focused on the illusory movement elicited by skin-stretch stimuli, when no actual motion was performed, often combined with local anesthesia to isolate cutaneous cues. Of note, these stimuli were typically delivered in a qualitative manner, leaving unclear how controlled, augmented skin deformation influences proprioceptive processing during active movement. In this study, we addressed this question by applying precisely scaled skin-stretch stimuli across the proximal interphalangeal (PIP) joint of the index finger, amplifying skin deformation patterns naturally occurring during finger flexion. Participants performed an active hand position-matching, with their bare hand and while wearing a custom-built, non-invasive tactile wearable device that can deliver controlled skin-stretch stimulation across PIP. While no significant differences in hand position were observed between the bare hand case and the condition where the device was worn but deactivated, augmenting the natural skin stretch consistently shifted the perceived finger posture toward more extended configurations. This finding demonstrates that cutaneous deformation can reshape proprioceptive estimates even when muscle spindles and efferent signals are fully engaged, revealing a continuous integration of tactile and muscular-skeletal cues in kinesthetic perception. These results advance our understanding of skin stretch contribution to proprioceptive processing in active tasks, opening interesting scenarios for haptics-engaged human machine interaction and extended reality. 

The experimental data we collected are publicly available here.


# Hand Matching Experiment
Twelve participants performed an hand-matching task designed to assess proprioceptive reproduction of finger posture under different sensory conditions. The task involved using the right hand to grasp a cylinder of specified diameter and replicating the posture with the left hand. The cylinders were 3D-printed in ABS with fixed diameters of 45, 60, 75, and 90 mm.
Each trial began with both hands open and aligned in the parasagittal plane, with the palms facing medially and the forearms resting on the table. The experimenter then placed the reference cylinder in the participant’s right hand, serving as the cue to close the hand and replicate the posture with the left hand. Once participants confirmed they had assumed the correct posture, they reopened both hands, and the experimenter proceeded with the next cylinder. Cylinder sizes were presented in randomized order across trials to minimize potential learning effects and ensure unbiased replication performance.
Participants performed the hand-matching task under three conditions. In condition A, the device remained inactive, In condition B, participants performed the task without wearing the device. In condition C it delivered movement-contingent skin stretch. In condition B, participants performed the task without wearing the device. PIP joint angles were recorded either via the encoder embedded in TWIST (conditions A and C) or using a Leap Motion sensor mounted on opaque goggles (condition B), which prevented visual feedback of the environment while capturing accurate joint kinematics. 

# Dataset Content

This dataset comprises measurements from 12 participants performing the hand matching experiment. Each participant is represented by a .mat file (sub1.mat … sub12.mat) containing a 3×4 cell array of PIP (proximal interphalangeal) joint angle trajectories. Each cell contains the average waveform across trials for the corresponding condition and cylinder diameter.

* **Rows (3)** - Experimental conditions:

    * *DevOff*: Device worn, motors off;

    * *BareHand*: No device;

    * *DevOn*: Device worn, motors on.


* **Columns (4)** - Cylinder diameters:

    * *D1*: 45mm;

    * *D2*: 60mm;

    * *D3*: 75mm;

    * *D4*: 90mm.


# Included Analysis 
The dataset includes LMManalysis.R, which performs statistical analysis of the plateau angle of the PIP joint. In brief, the analysis comprises:
1. Plateau extraction: For each waveform, the region surrounding the peak (plateau) where the PIP angle stabilizes is identified. The function extract_plateau() calculates the mean plateau angle, start and end indices, and plateau length.
2. Data aggregation: Plateau measurements are collected into a single tibble with columns: Subject, Condition, Cylinder, PlateauAngle, i_start, i_end, plateau_len.
3. Linear Mixed Model (LMM): A linear mixed-effects model lmer(PlateauAngle ~ Condition + Cylinder + (1|Subject)) evaluates the effects of experimental condition and cylinder diameter on PIP plateau angle, accounting for inter-subject variability as a random effect.
4. Post-hoc pairwise comparisons: Pairwise contrasts between conditions (DevOff vs BareHand, DevOff vs DevOn, BareHand vs DevOn) are computed using emmeans.
5. Visualization: Boxplots and bar plots illustrate the distribution of plateau angles across conditions and cylinder diameters.
Analyses and figures can be reproduced by running LMManalysis.R in R with the following packages: R.matlab, dplyr, lme4, emmeans, and ggplot2.



# Acknowledgment
This work was supported in part by the Italian Ministry of University and Research (MUR) - Fondo Italiano per la Scienza (FIS), with the grant PERCEIVING (no. FIS00001153); in part by the European Research Council Synergy Grant Natural BionicS (NBS) project (Grant Agreement No. 810346), and in part by the Next Generation EU Project ECS00000017 ‘Ecosistema dell’Innovazione’ Tuscany Health Ecosystem (THE, PNRR, Spoke 9: Robotics and Automation for Health).



