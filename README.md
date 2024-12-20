# Figure Ground Segregation in Texture Stimuli

MATLAB code utilized for data acquisition in the psychophysics experiment described in the article:


Karimian, A., Roberts, M.J., De Weerd, P., & Senden, M. (n.d.). Gamma Synchrony Mediates Figure Ground Perception. *Manuscript submitted*.

## Abstract
Gamma synchrony is ubiquitous in visual cortex, but whether it contributes to perceptual grouping remains contentious based on observations that gamma frequency is not consistent across stimulus features and that gamma synchrony depends on distances between image elements. These stimulus 
dependencies have been argued to render synchrony among neural assemblies encoding components of the same object difficult. Alternatively, these dependencies may shape synchrony in meaningful ways. Using the theory of weakly coupled oscillators (TWCO), we demonstrate that stimulus dependence is crucial for gamma's role in perception. Synchronization among coupled oscillators depends on frequency dissimilarity and coupling strength, which in early visual cortex relate to local feature dissimilarity and physical distance, respectively. We manipulated these factors in a texture segregation experiment wherein human observers identified the orientation of a figure defined by reduced contrast heterogeneity compared to the background. Human performance followed TWCO predictions both qualitatively and quantitatively, as formalized in a computational model. Moreover, we found that when enriched with a Hebbian learning rule, our model also predicted human learning effects. Increases in gamma synchrony due to perceptual learning predicted improvements in behavioral performance across sessions. This suggests that the stimulus- dependence of gamma synchrony is adaptable to the statistics of visual experiences, providing a viable neural grouping mechanism that can improve with visual experience. Together our results highlight the functional role of gamma synchrony in visual scene segmentation and provide a mechanistic explanation for its stimulus-dependent variability.

## Experiment Overview

Eight participants engaged in a task to discriminate the orientation of a rectangular figure within a texture stimulus made up of Gabor annuli. This figure region was distinguishable from the background by its contrast distribution, manipulated through two independent variables: **Contrast Heterogeneity** and **Grid Coarseness**.

- **Contrast Heterogeneity**: Refers to the variance in contrasts of the Gabor annuli within the figure region. This experiment utilized five levels of contrast heterogeneity for the figure: 0.01, 0.2575, 0.505, 0.7525, and 1. The background contrast was always set to maximum heterogeneity (1).
  
- **Grid Coarseness**: Represents the spacing between Gabor annuli, affecting the texture's density. The experiment tested five levels of grid coarseness: 1, 1.125, 1.250, 1.375, and 1.5.

The study comprised 9 sessions, each containing multiple blocks of trials. Every block presented 25 distinct combinations of Contrast Heterogeneity and Grid Coarseness.

### Stimuli, Tasks, and Procedure

- Each stimulus was a full-screen irregular grid of non-overlapping Gabor annuli on a grey background, with each annulus having a diameter of 0.7°, a spatial frequency of 5.7 cycles/degree, and a mean luminance of 60.76 Cd/m2.
- A rectangular figure differing in contrast distribution from its surrounding texture was embedded in the lower right quadrant for sessions 1-8 and moved to the upper left quadrant for session 9 (transfer session) to test learning transfer.
- Participants indicated the figure's orientation (horizontal or vertical) using arrow keys, within a two-alternative forced-choice setup. The stimulus remained visible for up to 1000 ms unless a response was made or fixation was lost.
- Correct responses were followed by a green fixation point (500 ms), incorrect by a red point.
- The experiment was conducted in a dimly lit room, using a chin and head-rest to maintain a constant eye-screen distance of 57 cm on a 19'' Samsung SyncMaster 940BF LCD monitor. Stimuli display and response recording were facilitated by Psychtoolbox-3 for MATLAB, with eye movements monitored by an Eyelink 1000 eye-tracker.

### Procedure for Handling Aborted Trials

Trials where participants broke fixation were aborted and rescheduled randomly within the session.

### Parametrization of the Experimental Setup

- Eye-screen distance: 57 cm
- Stimulus presentation time: 1000 ms
- Inter-trial interval: 900 ms
- Each session included 30 blocks, each with 25 trials.

### Transfer Session

Session 9 tested the transferability of learned skills by relocating the figure to the upper left quadrant, with participants uninformed of the specific quadrant to expect the figure.

### Ethics

Participants provided written consent after full briefing, in alignment with the Helsinki Declaration. The local Ethical Committee of the Faculty of Psychology and Neuroscience (ERCPN) approved all procedures. Participants were compensated for their involvement.

## Requirements

- MATLAB 64-Bit (Version 3.0.14 or later)
- Psychtoolbox-3
- Eyelink 1000 eye-tracker (for fixation monitoring)

## License

This project is open-sourced under the MIT License. See the `LICENSE` file for more details.
