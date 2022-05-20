.. _ci_jobs:

***************************************
What runs on our Continuous Integration
***************************************

Overall the aim of the Continuous Integration is to provide a way for contributors to test their code on multiple platforms and help in the development process.
While we can run many jobs, we aim to limit the selection depending in the context to ensure that the time taken is not wasted or ends up delaying a pull request or release.

The goal is that, the builds we do not run on a pull request should primarily be ones which are unlikely to be the fault of the change in the pull request if they fail.

Currently we have several stages of CI jobs, some run on a pull request and some will only run on a schedule.

1. "core" - Pull Request, Scheduled and Release
    This runs a basic offline test suite on Linux for the latest version of Python we support.
    It ensures that basic functionality works and that we don't have any regressions.

2. "test" - Pull Request, Scheduled and Release
    This runs a basic test suite on Windows and Mac OS for older versions of Python we support.
    It ensures that basic functionality works and that we don't have any regressions.
    This stage needs to wait for the "core" stage to complete.
    Furthermore, we run:

    * "oldestdeps" - Check the offline test suite with the oldest dependencies installed.

3. "docs" - Pull Request, Scheduled and Release
    This runs a documentation build (without executing gallery examples) to test that the HTML documentation can be generated without any errors or warnings.
    The build is cached and the gallery examples are then tested during the "online" stage.

4. "online" - Pull Request, Scheduled
    This runs a full test suite on Linux for a version of Python.
    This build can fail (due to a range of reasons) and can take up to 45 minutes to run.
    We want to ensure that the offline tests are passing before we consider running this.
    Therefore this stage needs to wait for the "core" stage to complete.

    In addition, we run the documentation build to execute the gallery examples which can fail due to the same reasons as the online build.
    The documentation build from the "docs" stage is cached and restored for this documentation build.
    Therefore this stage also needs to wait for the "docs" stage to complete.

    As these are not tested when we build the distribution or wheels, we skip these on a release.
    These should be checked based on the last commit for that release branch instead before a tag.

5. "cron" - Scheduled
    Here we put builds that are useful to run on a schedule but not so useful on a pull request.
    This allows us to run more "exotic" or focused builds that day to day, should not affect a pull request.

    These are:

    * "base_deps" - Check that sunpy does not error if you only have the base dependencies installed.
    * "devdeps" - Check the offline test suite with the development dependencies installed.
      Likely to break due to upstream changes we can't fix and have to wait to be resolved.
    * "conda" - Check the offline test suite when using conda instead of pip/pypi.
      Likely to break due to packaging upstream we can't fix and have to wait to be resolved.
