# Chain sequence for running scripts:
#   - FLKEYS_data_processing.R
#   - plots_obs.R
#   - output_basic.R
#   - plots_basic.R
#   - output_multi.R


#           Consider:

# Master script (e.g., master_script.R)

# Source the first script
source("script1.R")

# Source the second script
source("script2.R")

# Source the third script
source("script3.R")

# Continue with additional scripts as needed



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


# plotting consideration:
# - Make a plot showing all of the coral-level "SIR"s (but maybe just the infection curve? or a separate plot with the dead lines?), with the x-axis being days from 0 to max number of days a coral was infected for. this would maybe just be a nice way to visualize 
#       all of the distributions of infection, and show off the internal workings of the model and data preparation
# - And maybe a similar plot, but the x-axis being survey dates, so you can see visualize the rise and fall of every coral's infection curves throughout time. could plot everything and maybe color map by site
# - These could be nice accompanying figures to ones showing a range of simulations of predicted outbreaks


# considerations for code I need to dredge up:
#   - the SIR testing scripts - for multi, basic, tissue vs count, "human population", etc. there were a ton of these kinds of scripts at one point --> *this is in /archive*
#   - the code to run a range of simulations with my fitted values (or maybe other values ?) --> *at least some of this is in /archive*. this made nice plots. there is one plot version that may be missing; can visit Teams chat history from last December
#   - the code to test different curves for lambda --> *these might be on Box archive /archive/Chapter1_Epidemics/FLKEYS_SIR/misc/old_scripts*. Go through and systematically archive these properly!