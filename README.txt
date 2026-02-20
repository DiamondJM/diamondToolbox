================================================================================
                           sourceLocalizer — README
================================================================================

Diamond JM et al., "Travelling waves reveal a dynamic seizure source in human
focal epilepsy," Brain, 2021.

Diamond JM et al., "Interictal discharges in the human brain are travelling
waves arising from an epileptogenic source," Brain, 2023.

--------------------------------------------------------------------------------
OVERVIEW
--------------------------------------------------------------------------------

sourceLocalizer is a MATLAB class for localizing the cortical source of
epileptiform activity — either seizure (ictal) discharges or interictal
epileptiform discharges (IEDs) — from intracranial EEG (iEEG) recordings.

The core hypothesis is that both ictal and interictal discharges represent
travelling waves propagated outward from a relatively focal cortical source.
Rather than treating the electrode with the earliest discharge as the source
(a common clinical heuristic), this toolbox uses time or phase differences of
arrival across nearby electrode pairs to triangulate the source location on the
cortical surface. This approach is analogous to multilateration methods used in
acoustics, radar, and earthquake epicenter detection.

Localization does not require that any electrode be implanted directly at the
source — the source can be inferred from surrounding electrodes that receive
the travelling wave.

--------------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------------

- MATLAB (with the Signal Processing Toolbox for filtfilt, hilbert, findpeaks,
  buffer, medfilt1)
- diamondToolbox (this repository), added to the MATLAB path via
  addpath(genpath(pwd)) from the diamondToolbox root directory
- braindata2 and brainplotter objects (can be created with diamondToolbox), which
  handle cortical surface meshes and cortical surface reconstruction plotting
- A subject folder with the following expected structure:

      rootFolder/
          <subj>/
              tal/
                  leads.csv        — electrode localization file (see below)
              brainData_<subj>.mat — cached braindata objects (auto-created)
              geodesic_<subj>.mat  — cached geodesic distances (auto-created)

- leads.csv must contain columns: chanName, x, y, z
  (electrode name and RAS coordinates in mm)

--------------------------------------------------------------------------------
QUICK START
--------------------------------------------------------------------------------

1. Navigate to the diamondToolbox root directory in MATLAB.

2. Construct a sourceLocalizer object:

      sl = sourceLocalizer(subj, rootFolder, timeSeries, chanNames, Fs);

   Inputs:
      subj        — Subject identifier string, e.g. 'NIH032'
      rootFolder  — Path to the root data folder (see folder structure above)
      timeSeries  — [n x m] numeric array; n samples, m channels
      chanNames   — [m x 1] cell array of electrode name strings, 
                    where the elements of chanNames correspond to the columns of timeSeries
      Fs          — Sampling rate in Hz


   On construction, the object loads or creates braindata objects (myBd,
   myBp) for the subject, which contain the cortical surface mesh and
   associated utilities. These are cached to disk for reuse.

   Note: the chanNames column in leads.csv does not need to match chanNames. 
   leads.csv should contain every element of chanNames, but can additionally 
   contain extra elements. If there are elements of chanNames not present in 
   leads.csv, a warning will be thrown, and the missing channels cannot 
   participate in localization, but the code should proceed. 

   Critically, however, there should be consistency in channel naming between 
   chanNames and leads.csv. For instance, for some channel chanX, the time series 
   of chanX given in timeSeries must have been recorded at the location of chanX
   given in leads.csv. 

3. Set the localization mode (default is 'spikes'):

      sl.localizationMode = 'spikes';    % for interictal discharge localization
      sl.localizationMode = 'seizure';   % for ictal discharge localization

4. Run the full localization pipeline:

      sl = sl.localizationManager();

   This executes spike detection or seizure phase extraction, computes time or
   phase differences across sensor pairs, runs the geometric source
   localization, maps results to anatomical regions of interest, and produces
   summary plots.

5. Optional arguments for localizationManager:

      sl = sl.localizationManager('plotting', false);  % suppress plots
      sl = sl.localizationManager('forceNew', true);   % ignore cached results

--------------------------------------------------------------------------------
KEY PARAMETERS
--------------------------------------------------------------------------------

Default parameters are set automatically within prepareDeltaPosition and can
be modified after construction but before running localizationManager. They are
stored in sl.sourceLocalizationResults.paramStruct.

  propagationSpeed   — Assumed wave propagation speed in mm/s.
                       Default: 300 mm/s for both modes.

  sensorDistance     — Maximum geodesic distance (mm) between two electrodes
                       for them to be treated as a sensor pair.
                       Default: 30 mm (spikes), 12 mm (seizure).

  subsensorLength    — Minimum number of sensor pairs (or electrodes in a
                       sequence) required for a localization estimate to be
                       computed.
                       Default: 3 (spikes), 4 (seizure).

  distanceThresh     — Maximum geodesic distance (mm) from a candidate source
                       vertex to an electrode in the pair. Vertices beyond this
                       threshold are excluded from the candidate source space.
                       Default: 30 mm. Reflects the empirical limit of local
                       field potential spread (~3 cm) reported in the
                       literature.

For spike detection, parameters are stored in
sl.spikeDetectionResults.paramStruct:

  zThresh            — Z-score prominence threshold for peak detection.
                       Default: 3.
  ampScale           — The combined positive and negative peak prominence must
                       exceed ampScale * zThresh for a spike to be accepted.
                       Default: 3.
  maxNegPeakWidth    — Maximum width of the negative peak in seconds.
                       Default: 0.05 s.
  seqWin             — Window duration (seconds) used to group spikes into
                       sequences.
                       Default: 0.1 s.
  peakWin            — Window around a detected peak used to extract the
                       spike waveform.
                       Default: 0.1 s.

--------------------------------------------------------------------------------
PIPELINE DESCRIPTION — SPIKE (INTERICTAL) MODE
--------------------------------------------------------------------------------

Brain 2023 paper: IEDs are treated as travelling waves that propagate from a
focal epileptogenic source. Sequences of IED arrivals across electrodes are
used to estimate the source location.

Step 1: Spike Detection  (findSpikeTimes)
-----------------------------------------
Each channel of the input time series is z-scored. Negative peaks are
identified using MATLAB's findpeaks with a minimum prominence of zThresh
(default 3 SD) and a maximum width of maxNegPeakWidth. A spike is accepted
only when a nearby positive peak exists within peakWin seconds such that the
combined positive and negative prominence exceeds ampScale * zThresh.

Spikes are stored as a sparse binary raster (samples x channels).

Volume conduction removal: Any time sample at which two or more electrodes
simultaneously spike is flagged as likely volume-conducted and removed.

Step 2: Sequence Computation  (computeSequences)
-------------------------------------------------
The raster is buffered into overlapping 100ms windows (seqWin). Windows
containing spikes on at least subsensorLength distinct electrodes are
retained as candidate IED sequences.

Within each window, spike times are sorted by time of arrival. Duplicate
spikes on the same electrode within a window (e.g. from complex waveforms)
are resolved by selecting the time closest to the median.

Duplicate sequences — defined as sequences sharing a spike on the same
electrode at the same absolute time — are deduplicated, retaining the longest
version.

Step 3: Sensor Pair Construction  (getSensors)
-----------------------------------------------
Geodesic distances between all electrode pairs are computed on the pial
surface mesh (see Geodesic Distances below). Pairs whose geodesic distance
is within sensorDistance (default 30 mm) and on the same hemisphere are
retained as sensor pairs.

Step 4: Sequence to Delta Position  (spikeSequenceToDeltaPosition)
-------------------------------------------------------------------
For each IED sequence and each valid sensor pair, the difference in spike
arrival times (in samples) between the two electrodes in the pair is
recorded. This value, deltaPosition, is the raw observation from which
source location is inferred.

Only inter-electrode delays within the theoretically expected maximum
(sensorDistance / propagationSpeed) are retained.

--------------------------------------------------------------------------------
PIPELINE DESCRIPTION — SEIZURE (ICTAL) MODE
--------------------------------------------------------------------------------

Brain 2021 paper: Ictal rhythmic activity is treated as a travelling wave.
Instantaneous phase differences across electrode pairs are used to estimate
the source location continuously throughout the seizure.

Step 1: Phase and Power Extraction  (extractPhasePower)
--------------------------------------------------------
The time series is bandpass filtered using seizureWindowFilt (a fixed filter
targeting ictal frequencies, included in diamondToolbox). The Hilbert
transform is applied to extract instantaneous phase, amplitude, and frequency
(as the time derivative of unwrapped phase, smoothed with a 1000-sample
median filter).

Step 2: Phase Post-Processing  (postProcessPhase)
--------------------------------------------------
At each time sample, only the 14 channels with the highest instantaneous
amplitude are retained; the rest are set to NaN. This focuses the analysis
on electrodes with strong ictal signal and discards low-SNR channels.

Step 3: Phase to Delta Position  (phaseToDeltaPosition)
--------------------------------------------------------
For each sensor pair and each time sample, the phase difference between the
two electrodes is computed with explicit wrapping correction — the algorithm
selects the wrapped version of the second electrode's phase that minimizes
the absolute phase difference (i.e., it chooses between
phase2 - 2π, phase2, phase2 + 2π).

A Doppler correction is applied: pairs are excluded at time points where the
two electrodes' instantaneous frequencies differ by more than an acceptable
ratio, derived from the assumed propagation speed (300 mm/s) and an assumed
maximum source velocity (20 mm/s):

      acceptableRatio = (v + v_source) / (v - v_source)

The phase difference is then converted to a spatial offset (mm):

      deltaPosition = -(phi_1 / (2π * f_1) - phi_2 / (2π * f_2)) * v

where phi is the mean-subtracted phase at each electrode and f is
instantaneous frequency. Positive deltaPosition means electrode 1 is closer
to the source than electrode 2.

--------------------------------------------------------------------------------
CORE LOCALIZATION  (localizationFunction)
--------------------------------------------------------------------------------

This function is mode-agnostic; it operates on the deltaPosition matrix
produced by either pipeline above.

Geometric Framework
-------------------
For a given sensor pair (electrode i, electrode j), the measured
deltaPosition equals:

      geodesic_distance(source → i) − geodesic_distance(source → j) = ΔD

This defines a hyperbola on the cortical surface: the locus of all points
whose geodesic distance difference to the two electrodes equals ΔD. The
true source lies at the intersection of the hyperbolas from all sensor pairs
active at that time point.

Geodesic distances from every electrode to every vertex of the pial mesh are
precomputed and cached (see Geodesic Distances below).

Localization Procedure (per time point)
----------------------------------------
1. Active sensor pairs are identified — those with a non-NaN deltaPosition
   that is less than perceivedActualCutoff (0.9) times the inter-electrode
   geodesic distance. This cutoff excludes eccentrically-shaped hyperbolas
   that occur when the perceived delay approaches the full inter-electrode
   distance, as these provide poor geometric constraints.

2. Hemisphere assignment: the hemisphere used is determined by majority vote
   among the active sensor pairs.

3. Candidate source vertices are restricted to those within distanceThresh
   (30 mm) of at least one electrode in any active pair — the union of the
   sensor pairs' source spaces.

4. For each active sensor pair, the hyperbola is materialized as the set of
   candidate vertices satisfying:

      |geodesic(v → i) − geodesic(v → j) − ΔD| < marginError

   where marginError = 0.5 mm.

5. For each candidate vertex in the source space, the minimum Euclidean
   distance from that vertex to each hyperbola is computed.

6. The source is estimated as the candidate vertex that minimizes the mean
   squared distance (L2 norm) across all hyperbolas.

Quality Control
---------------
Estimates with a mean squared residual exceeding qualityControlThresh
(10 mm) are discarded. This removes localizations where the hyperbolas
do not intersect well, indicating inconsistent or noisy observations.

Output
------
Localization results are stored in:

      sl.sourceLocalizationResults.localizationResults

This is a [3 x K] matrix where each column is a retained localization:
      Row 1: vertex index (negative = left hemisphere, positive = right)
      Row 2: mean squared residual (localization quality, mm)
      Row 3: time index (sample number in the original time series)

--------------------------------------------------------------------------------
GEODESIC DISTANCES  (collectGeodesicDistances_master, loadGeodesic)
--------------------------------------------------------------------------------

Geodesic distances are computed on the pial surface mesh rather than in
Euclidean space. This is important because the cortical surface is highly
folded; straight-line distances between electrodes or between an electrode
and a candidate source vertex poorly reflect true anatomical proximity.

The procedure:

1. Electrode coordinates from leads.csv are projected onto the nearest vertex
   of the pial surface mesh (left or right hemisphere determined by x < 0).
   Electrodes more than 5 mm from the nearest pial vertex (e.g., depth
   electrodes that cannot be projected) receive a NaN vertex assignment and
   are excluded from localization.

2. For each successfully localized electrode, the geodesic distance from that
   electrode's pial vertex to all other vertices on the same hemisphere is
   computed using brainplotter's dist_geodesic method.

3. Contralateral distances are set to NaN, restricting all sensor pairs and
   source localization to a single hemisphere.

Geodesic distances are cached to disk as geodesic_<subj>.mat and reloaded
on subsequent calls.

--------------------------------------------------------------------------------
POST-LOCALIZATION: ROI MAPPING  (locDataToRoi)
--------------------------------------------------------------------------------

Localized vertex indices are mapped to anatomical regions of interest (ROIs)
using braindata2's vertex2ROI method. For each ROI that contains at least one
localized vertex, the number of localizations falling within that ROI is
counted.

Results are stored as a containers.Map object (vertexMap) in:

      sl.sourceLocalizationResults.roiResults

Left-hemisphere ROIs receive a negative key sign; right-hemisphere ROIs
receive a positive key sign.

An optional timeWindow parameter (in seconds) allows ROI counts to be
accumulated within sliding temporal windows rather than collapsed across the
full recording, enabling visualization of how the source evolves over time
(most relevant for seizure mode).

See Trotta, et al. Human Brain Mapping 2015. 

--------------------------------------------------------------------------------
VISUALIZATION
--------------------------------------------------------------------------------

plotDimensionsReducedWrapper  /  plotDimensionReduced
------------------------------------------------------
Displays localized source positions projected onto three 2D anatomical views
(axial, coronal, sagittal) within a single figure.

Each point represents a unique localized vertex. Point size is proportional
to local density — how many other localized vertices fall within 2 mm.
Points are sorted and plotted largest first so denser clusters are visible
beneath sparse ones.

Color coding (colorMode):
  'heatmap' (default for spikes)  — color encodes local density via the
                                    turbo colormap.
  'timing'  (default for seizure) — color encodes mean time within the
                                    recording at which that vertex was
                                    localized (parula colormap), allowing
                                    visualization of dynamic source movement.

The brain boundary is drawn from the combined left and right pial mesh
vertices, subsampled and projected to 2D.

plotSurfFun
-----------
Renders localized ROIs on a 3D pial surface brain, colored by count. If a
timeWindow was specified in locDataToRoi, this function animates the
localization as a video, frame by frame, showing how the source moves over
the course of the seizure.

--------------------------------------------------------------------------------
OUTPUTS SUMMARY
--------------------------------------------------------------------------------

After running localizationManager, the following fields are populated:

  sl.spikeDetectionResults
      .rasters          — sparse [samples x channels] binary IED raster
      .waveforms        — cell array of spike waveforms per channel
      .paramStruct      — spike detection parameters used

  sl.seqResults         (spike mode only)
      .seriesAll        — cell array of electrode sequences
      .timesAll         — numeric array of inter-electrode timing differences

  sl.seizureProcessingResults   (seizure mode only)
      .phase            — instantaneous phase [samples x channels]
      .pow              — instantaneous amplitude
      .freq             — instantaneous frequency

  sl.sensors
      .sensorInds       — [P x 2] array of sensor pair electrode indices
      .distancesMaster  — pairwise geodesic distance matrix

  sl.deltaPosition      — [P x T] matrix of inter-electrode spatial offsets
                          (mm after conversion); P = sensor pairs, T = time

  sl.sourceLocalizationResults
      .localizationResults — [3 x K] matrix of retained localizations
      .roiResults          — ROI-mapped localization counts
      .paramStruct         — all parameters used in this run

--------------------------------------------------------------------------------
NOTES AND CAVEATS
--------------------------------------------------------------------------------

- Hemisphere constraint: Localization is performed per hemisphere. Sensor
  pairs spanning hemispheres are excluded. The hemisphere used for each time
  point is chosen by majority vote among active sensor pairs.

- Depth electrodes: Electrodes that cannot be projected to within 5 mm of
  the pial surface (typically deep structures such as hippocampus) will
  receive NaN geodesic assignments and will not contribute to sensor pairs
  or localization.

- Propagation speed: The assumed speed of 300 mm/s is consistent with values
  reported in the literature for cortical travelling waves. This parameter
  can be adjusted, but the localization result is relatively robust to
  moderate changes in propagation speed because the hyperbolic geometry
  depends more on the difference in distances than on the absolute scale.

- marginError (0.5 mm): This is the tolerance used to determine whether a
  candidate vertex lies on a hyperbola. It should be checked against the
  vertex density of the pial mesh being used. If the mesh is coarse, this
  value may need to be increased to avoid excessively sparse hyperbolas.

- Volume conduction: The spike detector conservatively removes all
  simultaneous multi-electrode spikes. In recordings with high IED rates or
  spatially dense coverage, this may discard a significant fraction of events.
  The fraction removed is reported to the console.

- Parallel processing: The spike detector (findSpikeTimes) and the core
  localization loop (localizationFunction) use parfor. A MATLAB Parallel
  Computing Toolbox license and an active parallel pool will accelerate
  these steps substantially.

--------------------------------------------------------------------------------
CITATION
--------------------------------------------------------------------------------

If you use this toolbox, please cite:

Diamond JM, Diamond BE, Trotta MS, Dembny K, Inati SK, Zaghloul KA.
"Travelling waves reveal a dynamic seizure source in human focal epilepsy."
Brain. 2021;144(6):1751–1763. https://doi.org/10.1093/brain/awab089

Diamond JM, Withers CP, Chapeton JI, Rahman S, Inati SK, Zaghloul KA.
"Interictal discharges in the human brain are travelling waves arising from
an epileptogenic source." Brain. 2023;146(5):1903–1915.
https://doi.org/10.1093/brain/awad015

================================================================================