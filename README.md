
<h1 align="center"> Augmenting Natural Skin Stretch During Active Finger Motion Reshapes Hand Proprioception </h1>


# Abstract
There is strong evidence that skin stretch at joint level contributes to proprioception, i.e., the inner perceptual representation of body position and motion. Considering the dorsum of the hand and fingers, state-of-the art studies focused on the illusory movement elicited by skin-stretch stimuli, when no actual motion was performed, often combined with local anesthesia to isolate cutaneous cues. Of note, these stimuli were typically delivered in a qualitative manner, leaving unclear how controlled, augmented skin deformation influences proprioceptive processing during active movement. In this study, we addressed this question by applying precisely scaled skin-stretch stimuli across the proximal interphalangeal (PIP) joint of the index finger, amplifying skin deformation patterns naturally occurring during finger flexion. Participants performed an active hand position-matching, with their bare hand and while wearing a custom-built, non-invasive tactile wearable device that can deliver controlled skin-stretch stimulation across PIP. While no significant differences in hand position were observed between the bare hand case and the condition where the device was worn but deactivated, augmenting the natural skin stretch consistently shifted the perceived finger posture toward more extended configurations. This finding demonstrates that cutaneous deformation can reshape proprioceptive estimates even when muscle spindles and efferent signals are fully engaged, revealing a continuous integration of tactile and muscular-skeletal cues in kinesthetic perception. These results advance our understanding of skin stretch contribution to proprioceptive processing in active tasks, opening interesting scenarios for haptics-engaged human machine interaction and extended reality. 

The experimental data we collected are publicly available here.


# Hand Matching Experiment
Twelve participants performed a hand-matching task designed to assess proprioceptive reproduction of finger posture under different sensory conditions. The task required participants to grasp a reference object with the right hand and reproduce the same finger posture with the left hand. The reference objects were 3D-printed cylinders made of ABS plastic with fixed diameters of 45, 60, 75, and 90 mm.
Each trial began with both hands open and aligned in the parasagittal plane, palms facing medially, and forearms resting on a table. The experimenter placed the reference cylinder in the participant’s right hand. Participants then closed the right hand to grasp the object and actively matched the perceived finger posture with the left hand. Once the participant confirmed that the posture had been matched, both hands were reopened and the next trial began. Cylinder diameters were presented in randomized order to minimize learning effects.
Participants performed the hand-matching task under three conditions. In condition A, the device remained inactive, In condition B, participants performed the task without wearing the device. In condition C it delivered movement-contingent skin stretch. In condition B, participants performed the task without wearing the device. PIP joint angles were recorded either via the encoder embedded in TWIST (conditions A and C) or using a Leap Motion sensor mounted on opaque goggles (condition B), which prevented visual feedback of the environment while capturing accurate joint kinematics. 


# Supplementary Analysis
In addition to the main hand-matching experiment, supplementary data were collected to support additional analyses reported in the Supplementary Materials of the associated article. Specifically:
* Participants wore the device on the right hand while grasping the cylinders, allowing PIP joint angles to be recorded with the encoder and used as a reference baseline.
* During the BareHand hand-matching condition, the right hand was also tracked using the Leap Motion sensor while grasping the cylinders, providing a corresponding reference baseline measured with the same sensing modality.
These supplementary measurements enable direct comparisons between sensing systems and provide reference data for estimating accuracy and skin-stretch-related effects.


# Dataset Content
The dataset includes measurements from 12 participants performing the hand-matching experiment and the supplementary baseline measurements. Each participant is represented by a MATLAB file (sub1.mat … sub12.mat) containing a 5 × 4 cell array of PIP joint angle trajectories. Each cell stores the average joint-angle waveform across trials for a specific experimental condition and cylinder diameter.

* **Rows (3)** - Experimental conditions:

    * *DevOff*: Device on the left hand (matching), motors off;

    * *DevOn*: Device on the left hand (matching), motors on;

    * *BareHand*: No device, left hand trajectories;

    * *Baseline_Encoder*: Device on the right hand while grasping objects, motors off;

    * *Baseline_LeapMotion*: No device, right hand trajectories .


* **Columns (4)** - Cylinder diameters:

    * *D1*: 45mm;

    * *D2*: 60mm;

    * *D3*: 75mm;

    * *D4*: 90mm.


# Included Analysis 
This repository contains scripts that reproduce the main analyses reported in the article, as well as all supplementary analyses.

* **LMManalysis.R**: statistical analysis of the plateau angle of the PIP joint.
    1. Linear Mixed Model (LMM): A linear mixed-effects model lmer(PlateauAngle ~ Condition + Cylinder + (1|Subject)) evaluates the effects of experimental condition and cylinder diameter on PIP plateau angle, accounting for inter-subject variability as a random effect.
    2. Post-hoc pairwise comparisons: Pairwise contrasts between conditions (DevOff vs BareHand, DevOff vs DevOn, BareHand vs DevOn) are computed using emmeans.
    3. Visualization: Boxplots and bar plots illustrate the distribution of plateau angles across conditions and cylinder diameters.
    Analyses and figures can be reproduced by running LMManalysis.R in R with the following packages: R.matlab, dplyr, lme4, emmeans, and ggplot2.

* **1.1-EstimatedSkinStretch.R**: Estimation of average skin stretch (% strain) for the reference cylinder.
    1. Skin stretch estimation: Plateau angles are converted into estimated skin stretch values using a linear transformation based on predefined offset and gain parameters, with an additional gain applied in the DevOn condition, as defined in Methods Section.
    2. Normalization and summary metrics: Mean skin stretch is computed per condition and expressed as a percentage relative to a reference cylinder (Cylinder 1).
    All preprocessing and stretch estimation steps can be reproduced by running 1.1-EstimatedSkinStretch.R in R with the following packages: R.matlab, dplyr, tidyr, ggplot2, lme4, emmeans, and rstatix.

* **2.1-SkinStretchContribution.R**: Skin stretch contribution to the overall posture.
    1. Comparison between DevOff and DevOn: The mean plateau angles for reference cylinder are computed per condition, and the difference Δγ between DevOff and DevOn is calculated.
    2. Skin stretch estimation: Plateau angles are converted into estimated skin stretch using a linear transformation with condition-specific gain, and the difference Δa₀ between conditions is computed.
    3. Weight computation: The relative weight wC is calculated from Δγ and Δa₀ according to the derivative of the stretch–angle relation.
    All steps can be reproduced by running StretchEstimation_Cylinder2.R in R with the following packages: R.matlab, dplyr, tidyr, ggplot2.

* **3.1-SubjectVariability.R**: Intra-subject variability during grasping of real objects.
    1. Linear Mixed Model (LMM): A linear mixed-effects model lmer(PlateauAngle ~ Cylinder + (1|Subject)) evaluates the effect of cylinder diameter on PIP plateau angle, accounting for inter-subject variability as a random effect.
    2. Visualization: Boxplots with overlaid individual subject points illustrate the distribution of plateau angles across cylinder diameters.
    3. Summary statistics: Median, first quartile (Q1), third quartile (Q3), and interquartile range (IQR) are computed for each cylinder.
    All steps can be reproduced by running LMManalysis_Baseline.R in R with the following packages: R.matlab, dplyr, lme4, and ggplot2.

* **3.2-Bias.R**: Parametric and non-parametric comparison of baseline conditions.
    1. Primary analysis (parametric):
        * Linear Mixed Model (LMM): lmer(PlateauAngle ~ Condition + Cylinder + (1|Subject)) evaluates the effects of experimental condition and cylinder diameter on PIP plateau angle, accounting for inter-subject variability as a random effect.
        * Post-hoc pairwise comparisons: Contrasts between conditions are computed using emmeans.
        * Descriptive statistics: Median, first quartile (Q1), third quartile (Q3), and interquartile range (IQR) are computed separately for each condition and cylinder.
        * Visualization: Boxplots and bar plots illustrate plateau angle distributions across conditions and cylinders.
    2. Complementary analysis (non-parametric):
        * Friedman tests: Cylinder as a within-subject factor (per condition) and condition as a within-subject factor (per cylinder) are evaluated.
        * Post-hoc paired Wilcoxon tests: Pairwise comparisons with FDR correction are performed for cylinders within each condition and conditions within each cylinder.
        * Descriptive statistics: Median and interquartile range are computed for each condition × cylinder combination.
        * Visualization: Boxplots with subject-level points show distributions for the non-parametric analysis.
All steps can be reproduced by running LMManalysis_Comparative.R in R with the following packages: R.matlab, dplyr, lme4, emmeans, rstatix, and ggplot2.

* **4.1-Accuracy.R**: Comparison between hand-matching accuracy and augmented skin-stretch effects.
    1. Error computation:
        * Accuracy: difference between Baseline_LeapMotion and BareHand conditions, reflecting hand-matching error.
        * Augmented skin stretch effect: difference between Baseline_Encoder and DevOn conditions.
        * Residual difference: difference between skin stretch effect and accuracy.
    2. Descriptive statistics per subject:
        * Mean and standard deviation of accuracy.
        * Mean and standard deviation of skin stretch effect.
    3. Statistical analysis:
        * Normality check: Shapiro–Wilk test per cylinder on the residual difference.
        * Paired t-tests: compare skin stretch effect vs accuracy per cylinder, reporting mean ± SD, t-value, degrees of freedom, and p-value.
    4. Visualization: Bar plots of mean ± SD per cylinder, comparing accuracy and skin stretch effect, with clear labels and color coding.


# Acknowledgment
This work was supported in part by the Italian Ministry of University and Research (MUR) - Fondo Italiano per la Scienza (FIS), with the grant PERCEIVING (no. FIS00001153); in part by the European Research Council Synergy Grant Natural BionicS (NBS) project (Grant Agreement No. 810346), and in part by the Next Generation EU Project ECS00000017 ‘Ecosistema dell’Innovazione’ Tuscany Health Ecosystem (THE, PNRR, Spoke 9: Robotics and Automation for Health).



