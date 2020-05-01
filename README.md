
# bh-tsne

Barnes-Hut t-SNE

**This code was derived from** [T-SNE-Java](https://github.com/lejon/T-SNE-Java)

See [LICENSE.md](LICENSE.md) for further information!

## Maven coordinates:

**TODO** The library will soon be available at 

    <dependency>
        <groupId>de.javagl</groupId>
        <artifactId>bh-tsne</artifactId>
        <version>0.0.1-SNAPSHOT</version>
    </dependency>


## Project structure

The project consists of two modules:

- `bh-tsne`: The actual library offering the t-SNE functionality. It has
  the `ejml` (Efficient Java Matrix Library) as its only dependency
  
- `bh-tsne-demo`: A small application that allows running the t-SNE on
  different test data sets, and visualizing the results: 

![bh-tsne-screenshot-001.png](/screenshots/bh-tsne-screenshot-001.png)
  

## Changes done to the original code:

- Replaced the use of `ThreadLocalRandom` with a `Random` instance
  that is initialized with a fixed (but modifiable) random seed. 
  Since t-SNE is often used for scientific publications, having the 
  possibility to create the same result repeatedly is *crucial* for 
  reproducability. I also used this to verify (empirically) that all 
  subsequent changes did not affect the validity of the results.
  
- Collected all classes that are required for the `BHTSne` into a 
  single package - and *only* these classes!
  
- Made `private` what could be `private`. Made *default-visible* 
  what could be *default-visible*. Made `static` what could 
  be `static`. Removed all code that was then no longer used
  (which was most of it) 

- Removed unnecessary dependencies (for example, JAMA was used for a 
  single line of code that *was not called at all*). Specifcially:
  Removed *all* dependencies except for the one to EJML, which 
  has been replaced with the proper one in its latest available
  version.
 
- Replaced the `System.out.print...` calls by some (pragmatic)
  logging and progress report
 
- Added a trivial parallelization that brought roughly 40% speedup by 
  changing a few lines of code (and seems to be faster than the oddly 
  complicated `ParallelBHTsne` of the original implementation...)
  
- Allowed the computation to be interrupted by calling `interrupt()` on
  the executing thread
  
- Significantly reduced the number of memory allocations, mainly by 
  introducing the `DoubleArray` interface

- Offered the whole functionality in a single public class, with some
  basic JavaDoc
  
- TODO Offered the library in Maven Central
  

  