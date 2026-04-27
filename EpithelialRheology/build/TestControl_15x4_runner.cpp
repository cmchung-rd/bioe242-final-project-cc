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
#include "/workspace/EpithelialRheology/test/TestControl_15x4.hpp"

static TestControl_15x4 suite_TestControl_15x4;

static CxxTest::List Tests_TestControl_15x4 = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestControl_15x4( "/workspace/EpithelialRheology/test/TestControl_15x4.hpp", 131, "TestControl_15x4", suite_TestControl_15x4, Tests_TestControl_15x4 );

static class TestDescription_TestControl_15x4_TestControl_rep1 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestControl_15x4_TestControl_rep1() : CxxTest::RealTestDescription( Tests_TestControl_15x4, suiteDescription_TestControl_15x4, 370, "TestControl_rep1" ) {}
 void runTest() { suite_TestControl_15x4.TestControl_rep1(); }
} testDescription_TestControl_15x4_TestControl_rep1;

static class TestDescription_TestControl_15x4_TestControl_rep2 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestControl_15x4_TestControl_rep2() : CxxTest::RealTestDescription( Tests_TestControl_15x4, suiteDescription_TestControl_15x4, 371, "TestControl_rep2" ) {}
 void runTest() { suite_TestControl_15x4.TestControl_rep2(); }
} testDescription_TestControl_15x4_TestControl_rep2;

static class TestDescription_TestControl_15x4_TestControl_rep3 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestControl_15x4_TestControl_rep3() : CxxTest::RealTestDescription( Tests_TestControl_15x4, suiteDescription_TestControl_15x4, 372, "TestControl_rep3" ) {}
 void runTest() { suite_TestControl_15x4.TestControl_rep3(); }
} testDescription_TestControl_15x4_TestControl_rep3;

#include <cxxtest/Root.cpp>
