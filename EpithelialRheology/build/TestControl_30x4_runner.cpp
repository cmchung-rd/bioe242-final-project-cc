/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "/workspace/EpithelialRheology/test/TestControl_30x4.hpp"

static TestControl_30x4 suite_TestControl_30x4;

static CxxTest::List Tests_TestControl_30x4 = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestControl_30x4( "/workspace/EpithelialRheology/test/TestControl_30x4.hpp", 70, "TestControl_30x4", suite_TestControl_30x4, Tests_TestControl_30x4 );

static class TestDescription_TestControl_30x4_TestEpithelialGrowthInRectangle : public CxxTest::RealTestDescription {
public:
 TestDescription_TestControl_30x4_TestEpithelialGrowthInRectangle() : CxxTest::RealTestDescription( Tests_TestControl_30x4, suiteDescription_TestControl_30x4, 74, "TestEpithelialGrowthInRectangle" ) {}
 void runTest() { suite_TestControl_30x4.TestEpithelialGrowthInRectangle(); }
} testDescription_TestControl_30x4_TestEpithelialGrowthInRectangle;

#include <cxxtest/Root.cpp>
