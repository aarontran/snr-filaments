Meeting Thursday 2014 July 3
============================

A spontaneous drop-in with Brian, around 11a-12p.  Covered a good amount of
stuff, clarified some understanding of the physics of Ressler's paper, and
added some things to the to-do list (esp. for documentation/writeup).

After reading Sean's email, I put together a few thoughts...


## Obtaining `m_E` from analytic expression, using best fit values

Basically, using Sean's codes, fit energy/widths and get best fit `B_0` and
`\eta` values.  Then we could punch the results into equations (13), (23).
See my email to Sean -- my idea is that these fits should incorporate all the
information about energy-width scaling.

Brian: How do you know what mu is, though?  Because you have a whole bunch of
values for that.  It's worth checking/trying.  The values had better fit
nicely in a least squares sense, or else you're in trouble!


## Interpretation of `m_E`, `mu`, other outputs

if we have information about eta, B0, `m_E`, can we say something about which
values of mu are best, or give most reasonable results, or something???

Brian: how did they constrain mu?  Good question...

Brian had some questions regarding diffusion coefficient / other things to
clarify... but they might have been encompassed by my questions for Sean.


## Tycho shock speeds

Shock velocities in Tycho vary by up to a factor of 2x (Williams et al.,
2013, ApJ).  Refer to Brian's paper, and use the appropriate values of the
shock speed (and shock speed / 4 for the speed downstream of the shock).
Angles are measured starting from north, going east (counter-clockwise looking
at the sky).  Basically take your regions, and find the closest shock velocity,
interpolating as needed.


## Mundane things

Yes, I'm not insane, equations (13) and (23) do give positive results (probably
just a sign error).

Water boiler and printer?  ... whatever (ask Pam?)


Agenda
------

* Send Brian Sean's code/guide (done)
* Incorporate Tycho's variable shock speed into results
* Eventually, generate a flowchart of what calculations were done, what
  equations were used, etc... would have been nice to do similar thing for
  Ressler paper, really.  What we're retracing now.
* See Brian's email: generate a handout of region pictures (zoomed in),
  spectra with fits and numbers, profile fit figures.  Need to send to Sean and
  Steve so they can see what we're doing.

