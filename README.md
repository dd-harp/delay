# delay

## notes
  * in the SEIR with delayed E->I reaction there are actually only 2 "events", as once S->E occurs, then E->I is guaranteed to occur. S->E occurs with the standard rate. But once the person is deposited in E, they stay there for a deterministic time, tau.
    * The punchline is that once S->E fires, calculated in the standard way, we assign t + tau as the time the reaction completes. We do it in 2 shots.
  * This gets awfully confusing if there are competing risks (aka; a mu flow of people out of each compartment from death.)
    * In that case the "purely delay approach" of Barbuti et al. may be what we want. The problem is that its not very elegant and they don't explain the math behind it. If they even know it. I like Anderson and Thanh because they explain the process that is being sampled.
    * It may be the case that if they can die while in E we can still use Anderson. Update the state in two shots as before.
    * The thesis "Stochastic Population Dynamics with Delay Reactions" may have the answer
