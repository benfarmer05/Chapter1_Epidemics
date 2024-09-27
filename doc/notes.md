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
  
