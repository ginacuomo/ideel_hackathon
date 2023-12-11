## File/Data Description

The `data_out` directory contains modelled outputs for external partners at the first administrative unit. 

### Modelled risk outputs

The first set of outputs provide the estimated time that artemisinin resistance (ArtR) will 
increase from 5% to 25%, with assumed starting partner drug reistance equal to 5%. Where times 
are predicted to take longer than 60 years, these have been censored (">60"):

1. `central_times.csv`
1. `optimistic_times.csv`
1. `pessimistic_times.csv`

The second set of outputs provide the estimated time from now (based on spatial maps of current reistance prevalence)
for ArtR to reach 25%.

Similarly to before, where times are predicted to take longer than 60 years, these have been censored (">60"):

1. `central_times_prospective.csv`  (STILL TO BE COMPLETED)
1. `optimistic_times_prospective.csv`  (STILL TO BE COMPLETED)
1. `pessimistic_times_prospective.csv`  (STILL TO BE COMPLETED)

In all the outputs above, the **central** times provide the estimated times based on the central parameter estimate for each of the parameters that we explored and that are known to impact the speed of selection of ArtR (malaria prevalence, treatment related parameters, current resistance prevelance, types of ACT in use). 

The **optimistic** times assume the upper or lower estimate value for each parameter (depending on the direction of its effect on selection) such that the selection of ArtR will increase at its slowest. Conversely, the **pessimistic** times assume the upper or lower value for each parameter (depending on the direction of its effect on selection) such that the selection of ArtR will increase at its fastest.

Lastly, the uncertainty provided in each of the files, reflects the uncertainty that arises from our use of a stochastic, individual based model for the selection of ArtR. 
