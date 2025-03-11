# Chain sequence for running scripts:
#   - FLKEYS_data_processing.R
#   - plots_obs.R
#   - output_basic.R
#   - plots_basic.R
#   - output_multi.R
#   - plots_multi.R
#   - error_eval.R
#   - tables_figures.R


#           Consider:

# Main script (e.g., main_script.R)

# Source the first script
source("script1.R")

# Source the second script
source("script2.R")

# Source the third script
source("script3.R")

# Continue with additional scripts as needed



- from processing script:
  #   - now storing patient zero tissue and counts. this should help with understanding which initial conditions allow a model to "work"
  #   - would still be ideal to zero in on a standard initial condition, but may just not be possible with how SCTLD functions. might
  #       require using cover and N.site as scalars to define patient zero starting values

27 Sep 2024

- where did the code I used to create a range of simulated output from different coral cover scenarios go? would be very helpful to track down that script. can check chatgpt conversations, Box history, etc.

- from output_basic_COVER script:
-	  # some thoughts:
  #   - this seemed to work better when I was working with numbers really close to, but below, 1.0. struggling
  #       when working with values close to, but ABOVE, 1.0.
  #   - why might that be?
  #   - I was getting R0's that made more sense on January 26th, and also earlier in the afternoon Jan 29 (check Box version control)
  #
  #   - I tested something below where 1.0 is the limit for a coral cover of 100%...but not the same as what worked well before.
  #       - with a lambda of 1.0, you end up with a very large effect of coral cover on transmission. this means that wonky things can
  #           happen
  #   - theoretically, should there not be some very large lambda value for which coral cover explains so much of transmission that
  #       a cross-site transfer of beta from Nearshore to Offshore / Midchannel simply results in zero outbreak?
  #         - that's what I was seeing before (like with the sqrt function directly), where past certain thresholds, there is either
  #           zero outbreak after the transfer, or an outbreak that reduces to exactly the Nearshore outbreak, proportionally
  #             - I should probably get this back to that place. and maybe hit a stopping point and move on to the multi-group model again
  #   - okay! I sort of hit that point again: a lambda of 3.0 does indeed result in no outbreak for Midchannel & Offshore
  #
  #   - one takeaway I've got is, it seems like it might be really important to have values close to 1.0 so that the infection
  #       can take off quickly early in the outbreak. the total removal can get close with current methods below when transferring
  #       between sites, but there's usually a phase difference (temporal mismatch). consider looking at the other functions again, like
  #       a simple square root or the exponential one Dan did. these might provide clues too
  #         - part of the problem is that when I change lambda, it's changing the relationship of ALL sites to cover, but also
  #             their relationship with each other
  
    #     I think it would be smart to start migrating everything over to the style set up in 'temp.R'. this should make it easier to
    #     build off of main list of functions from which to pull and edit for basic SIR vs grouped SIR. that was always an issue before
    # - and then it's a matter of making sure the SIR output still makes sense! assuming it does, make sure it looks good for the grouped
    #     SIR as well
    # - from there, need to bring in effect of DHWs. consider writing 1-2 hr/day M/T/W to get the bones together, and then really work
    #     hard on Thursday to get a narrative together for the discussion section
    #
    # - other general thoughts:
    #     - would be nice to return to the plots of varying coral covers and compositions. confused about where that went
    #     - try to fit all three sites with one algorithm?
    #     - more than anything, focus on simplicity and what will get a paper out. can always add more analyses down the road, but it is
    #         critical to get this submitted to Ecology ASAP




- from plots_basic_COVER script:
  #to-do:
  # - optimize the prediction of nearshore parameters onto offshore. and vice versa probably
  # - keep testing effect of simply not fitting to infections at all in the first optimization step
  # - also, thinking about the Dobbelaere approach. I could try a fit that accounts for all site.loops summed together (or averaged)
  #     - what would the prediction look like? similar to the 2020 paper? there are important differences in how they accounted for area
  #     - could allow beta to do whatever it wants, and then calculate a gamma that a constant R0 dictates...does this make sense?
  #         - might allow the actual shape of the curve to change between site.loops
  #         - could also always test that more with the multi-scale models or a UVI-specific wheels study...eventually
  #         - also though, my theory is really that only transmission (beta) itself should change with site.loop density. so I like what I did, too
  #
  #     - our model is tissue-based too - which is a huge plus. while its predictions aren't perfect, it's probably a much closer reality
  #         than just assuming a bunch of whole coral colonies die
  #         - and ours is multi-species!
  #   - last thought: why not throw together a simple statistical model relating initial cover (& SA) with final remaining tissue? this
  #         could be a pretty interesting natural extension of the work, and prove whether what we set out to do could simply be replicated
  #         with something like a GAM
  
    - on that point, there are some good sources to consider (got these from chatGPT):
        - Miller, J. et al. (2009). "Coral disease: A potential link between ecological degradation and disease susceptibility." Disease of Aquatic Organisms, 87(2), 119-128.
        - Eakin, C. M., et al. (2010). "Coral bleaching: Patterns, processes, and causes." Coral Reefs, 29(2), 307-317.
        - De La Torre, C. et al. (2018). "Bayesian hierarchical modeling for disease dynamics in coral reefs." PLoS ONE, 13(2), e0192808.
        -Giri, K. et al. (2021). "Machine learning models for assessing coral health using remote sensing data." Marine Ecology Progress Series, 663, 139-153.
        - Ainsworth, T. D., et al. (2016). "Climate change disables coral bleaching protection on the Great Barrier Reef." Nature Climate Change, 6(1), 82-87.
        
        
        
        
# notes for paper:
#   - I accidentally was starting the epidemic at Nearshore with a polyp_SA associated with midchannel instead - and I think it caused the epidemic to start with too little tissue. and it's really interesting, the model will still "fit" but it takes way too long for infections to start really kicking off exponentially, and also too long for removal to kick in. a reminder that classic SIR models are very very sensitive to the amount of starting infection when used with a system like coral tissue
#   - The way I am predicting outbreaks from initial conditions, the *shape* of the outbreak is not dynamic. this probably is okay not to worry about when publishing, but it's probably an area that future work should zero in on. because it seems like the shape of the
#       outbreak really *should* be changing site to site
#           - Sort of along the lines of the "simple" statistical model I proposed above. my method may get us fairly close to a crude solution for predicted tissue lost to SCTLD - but is it really helping us understand the outbreak dynamics? I would say only a little
#   - The DHW data provides a unique opportunity to run some simple ancillary stats that test the effect of DHWs / SST on aspects of individuals corals' infection trajectories. not sure if this is worth it or not
#   - It will be nice to go back to the 'cover.csv' data and compare model-predicted loss of cover with that which was observed. we could make the assumption that no SCTLD-unaffected coral tissue was lost, and add that SA back into the output data before a final
#           conversion with the 1:3 ratio
#   - worth remembering that while susceptible tissue can only go down in my model, susceptible hosts as a whole can rise - may be worth plotting this. even though those hosts didn't really "recover" and theoretically any host is still susceptible even as it is
#           actively infected, because of the nature of multifocal lesions

# plotting consideration:
# - Make a plot showing all of the coral-level "SIR"s (but maybe just the infection curve? or a separate plot with the dead lines?), with the x-axis being days from 0 to max number of days a coral was infected for. this would maybe just be a nice way to visualize 
#       all of the distributions of infection, and show off the internal workings of the model and data preparation
# - And maybe a similar plot, but the x-axis being survey dates, so you can see visualize the rise and fall of every coral's infection curves throughout time. could plot everything and maybe color map by site
# - These could be nice accompanying figures to ones showing a range of simulations of predicted outbreaks


# considerations for code I need to dredge up:
#   - the SIR testing scripts - for multi, basic, tissue vs count, "human population", etc. there were a ton of these kinds of scripts at one point --> *this is in /archive*
#   - the code to run a range of simulations with my fitted values (or maybe other values ?) --> *at least some of this is in /archive*. this made nice plots. there is one plot version that may be missing; can visit Teams chat history from last December
#   - the code to test different curves for lambda --> *these might be on Box archive /archive/Chapter1_Epidemics/FLKEYS_SIR/misc/old_scripts*. Go through and systematically archive these properly!


# Dan discussion notes 23 Oct 2024:
#       - try letting the scalars fit when interacting with JUST sst (no threshold etc.). but let the baseline be the temperatures around when      #           SCTLD first hit Offshore site (like 28 / 29C). below that baseline, transmission can increase (non-linearly? who knows but try it),
#           and above that baseline, transmission crawls slower. can even try some extreme tests like turning beta and/or gamma to zero
#           under certain conditions, but for now just see if this approach alone can "improve" on what we have. I have a feeling it will
#           do the same thing I've already tried honestly, but worth a shot
#       - can still also test DHW - this would directly address the idea that bleaching and symbiont expulsion is involved. DHW = 4 is #     #             bleaching, and may be an interesting value to test
#        - think about observation-informed modeling - we know that SCTLD could start at Offshore at 29C, so assume that is our baseline.
#               kind of arbitrary, but it shows that we acknowledge SCTLD is capable of being active at that temperature OR higher
#        - Dan had some concerns about fitting with all the data available, including the "second wave" - that said, I think this approach makes sense because it still fits about the same as not doing so (but check this again), and that further shows that given a null hypothesis of no effect of temperature, the SIR model will produce a first wave that roughly matches the observations - can go from there to explicitly build in the effect of temperature and bleaching


# Notes 8 Nov 2024:
- Getting closer, I think, with bringing in a second wave. But, noticing that I should not be fitting the final timepoint of infections for *any* simulation. That is likely introducing problems, and may even be part of why the "first wave" looks so clean in the current manuscript. The model may be forcing the simulated infecteds down to zero too early! Try a couple things: 
1.) run analysis as-is, but dropping the last timepoint for infecteds
2.) run analysis for *just* first wave, dropping everything after 30.5C (and/or DHWs exceeds ~2.0 around 1 July 2019). see how this affects all my rates
3.) after the above steps, again introduce effects of SST / DHW for a "second wave" and see what happens
4.) ALSO, just noticed I am only fitting to removal for the basic SIR model!! need to address that. although it's doing remarkably well as-is actually
 # - the chatGPT thread is pretty helpful from today


# Notes 3 Feb 2025:
- Figures are coming together; likely just 3-4 total figures. One if them is non-code based (Figure 1; conceptual SIR diagram)
- Figure 2 is observations of infectious tissue. This is mostly finalized
       - but need to rework the top row to be proportional to 
 - Figure 3 is a composite figure demonstrating the fit, and projected fit, of Nearshore & Offshore. Column 1 is single-host model; Column 2 is multi-host model. May have column 3 for single-host model with thermal stress but not sure ?
        - need to double-check where the fit is or isn't getting cut off when thermal stress threshold is breached, and plot that in with an abline
       - may consider doing some proportional plotting here, too. particularly for column 2. the visuals might be identical to as-is, but the y-axis would be more informative and demonstrate the proportionality of each susceptibility group to the total tissue demonstrated by column 1 (max possible y-value is 1.0, but would be lower and match the top of the curves)
       - figure out what to do with the thermal stress before/after. how best to visualize?

# Notes 6 Feb 2025:
- Meeting today with Dan; decided:
- Include 'wave' column for error_eval table which tracks R-squared not just for the where the fit occurred (pre- or post-thermal threshold breach), but also for the entire outbreak (compare not just apples to apples, but also apples to oranges)
- fig 3; maybe make a solid line turn to dotted after thermal threshold is breached / fitting stops or starts. and/or vertical line / shaded region
    - also maybe symbols become red or open afterwards? symbols should also be smaller and perhaps have some transparency
- finally, make a projected single-host DHW nearshore to offshore and just see what happens !

# Notes 11 March 2025:
  #   - so ... last steps really are:
  #       1.) making sure if I want to stick with R-squared or absolute residual optimization
  #       2.) whichever I choose for #1, make sure it looks right in the output
  #       3.) triple check that calculation of R-squared and optimization is being done correctly
  #       4.) choose how to present effect of thermal stress (within existing figures, also maybe separate fine detail of it into fig. 4)
  #       5.) as appropriate, embed R-squared values into figures
  #       6.) work with Gaby to wrap up fig. 1
  #       7.) work with Dan to finalize typesetting on figs. 2-3 (or 2-4)
  #       8.) update manuscript accordingly to all above. ensure any necessary references to SIR methods are included
  #       9.) better plot out and integrate thermal stress
  #       10.) make top row of fig. 2 proportional
  #                 - - "proportion of infected tissue SA out of total susceptible SA", with a free y-axis, is useful because it will essentially show what we are showing already in           the top row, but the y-axis will be more informative because it will provide a peek at the prevalence. And that's important since Nearshore not only had a lot more infected SA than the other sites, but also much higher prevalence
- That goes back to one of the central messages of the paper, that coral density really matters. Nearshore may have blown up largely because of how tightly packed in its coral tissue SA was, relative to the surveyed reefs' size
  #       11.) fig. 2 should have x-axes in both rows
  #       12.) For eval, consider within- and across-group summation of RMSE
  #       13.) Maybe Table (supp. ?) with parameter definitions, table with R2/RMSE, table with boundary conditions and optimized parameter (and R0!) results?
  #       14.) update manuscript to reflect fig 1 changes to equations
  #       15.) note in paper how multi-host convergence is not achieved ... more so a range of outputs that all tend to "look" around the same in the end. something to explore in my chapter 3
  #        16.) note in paper the starting tissue value and its possible importance!
  #        17.) try higher itermax with single-host model ... try and figure out if there is a better lambda for it and multi-host too. though I think iterations won't matter for single-host really at all. converges to a consistent solution near every time. and actually also I just figured out I'm pretty sure I manually found the exact right lambda ... wow (1.0). so, now it's a question of why the multi-host projection performs better than the single-host one I think. a really cool result there. but also, how good is the "curve" I'm using to set lambda (and in turn modify the transmission rate) in the first place?

