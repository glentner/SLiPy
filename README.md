Astronomy Python Library
========================

This is a Python package containing modules I've developed to speed 
up my work flow. It mostly just contains functions for doing spectroscopy. 
It will continue to be developed. I'm sharing it in the hope that others 
may find it helpful. 

**Dependencies:** Python 3.x, astropy, matplotlib, numpy, scipy

Quick note: the subpackage **astrolibpy** was not developed
by me. It was coded by Sergey Koposov (@segasai) at Cambridge (then at least). 
I found it useful for performing velocity corrections on my spectroscopic 
data. I've modified several modules such that it can be imported and used in 
Python 3.x. See his README file.

## Contents:

* [**Astro:**](##Astro/) 
A subpackage containing modules pertaining to astronomy and 
spectroscopy. Modules include **Fits**, **Simbad**, **Plot**, 
**DataType**, **Correlate**, **Telluric**, and **Velocity**.

* [**Framework:**](##Framework/)
A subpackage containing general framework code that I use in the packeage. 
Modules include **Arguments**, **Options**, **Command**, and **Display**.

* [**astrolibpy:**](##astro/)
As stated previously, this package was not developed by me, only modified. I'm
not going to document it's usage here. Its name is unfortunate for me as it
is a bit over done with the convention I was already using, but for consistency
I will keep it as it was from the author.
