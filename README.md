# interference-cancellation
Snap Blind and Semiblind Multiuser Detection and Interference Cancellation

Interference cancellation (IC) is the strategy for forming an estimate of incurred interference, like intersymbol interference (ISI), co-channel interference (CCI), adjacent channel
interference (ACI), etc., and subtracting it from received signals before detection. Though it has intensively been investigated over the last decades, its commercial applications
for mitigating multiuser access interference (MAI) in mobile network can only be found in the last serval years [1]. Several IC schemes, such as successive IC and subspacebased IC, have been reported to be implemented for 3G
mobile network like cdma2000 and WCDMA [2]. The recent standardizing Femto cell further brought more interests and potential applications of IC techniques [3]. Compared with
other detection strategies, IC focuses more on interference estimation. And different interference estimation methods may lead to different IC schemes [4], [5], e.g. successive IC,
multistage detection, decision-feedback IC (DFIC) [6], [7], etc. DFIC, including minimum mean squared error (MMSE) DFIC [6] and decorrelating DFIC [7], belongs to the decisiondriven schemes that combines features of successive IC and
multistage detection [4]. Conventional IC receivers are known to be able to solve the near-far problem with the knowledge of the signature information of all users [4]. However this
assumption isn’t consistent with many practical situations where the receiver may only know the signatures of the expected signals not interfering signals. Recent research has
been devoted to semiblind/blind implementation of IC for the practical applications in which adaptive filter techniques, e.g., Wiener filtering [8], Kalman filtering [9] and subspace-based
implementations [10], are among the popular choices.


## Conventional Multiuser Access Signal Model
In the procedure of advanced multiuser receiver development, it is known that a proper received signal model can help us understand received signals as well
as receiver design. There are two popular multiuser signal models which have been intensively discussed for multiuser receiver design. They are the conventional multiuser signal model and the subspace-based multiuser
signal model. In the conventional signal model, each received signal is directly taken as a linear combination of actual signal signatures [1], [3], [6]. Most related blind multiuser receivers are developed either by explicitly estimating the signal signature [4] or by removing
interfering signal components using adaptive filtering techniques, e.g., the blind receiver design with Wiener filter [3] and Kalman filter [6] techniques. Though the
conventional signal model provides us a natural view of received signals, the involved signature waveforms or amplitudes information is unknown and it usually take
the receiver lots of efforts to obtain it before detection. 

## Subspace Multiuser Access Signal Model
For compensating the weakness of the conventional signal model, the subspace signal model is proposed with subspace-based signal processing techniques [5]. In the subspace signal model, each received signal is
taken as a linear combination of signal subspace bases, which can be obtained by subspace signal processing techniques on the autocorrelation matrix of received signals. Subspace signal mode can be taken as a result
of parametric signal modelling and provides a in-depth comprehension of received signals. Though subspacebased approaches don’t need explicitly estimate each user’s signal signature and the initialization and adaptive
speed are improved with good performance, the signal subspace formation procedure still is not trivial.

## The Proposed Blind Multiuser Access Signal Model
It is known that the conventional signal model provides us the foundation for both optimal and conventional multiuser receiver design and subspace signal model helps us understand signal underneath structure.
However, neither of them is easy enough for developing the blind multiuser receivers for high-speed CDMA systems [2]. In order to solve the near-far problem with
minimum prior knowledge and computation complexity, we propose a new blind multiuser model with directly connecting the current received signal and several previous received signal while no explicitly signal structure
estimation. With this blind signal model and widely employed signal estimation criteria including least squares (LS), minimum mean-squared errors (MMSE) and maximum likelihood (ML), several novel blind multiuser
receivers are developed. There is no statistical signal estimation or subspace separation procedure required. Only a minimum number of previously received signals and the desired user’s signal signature waveform and
timing are required. Hence the computation complexity and detection delay can be much reduced.
