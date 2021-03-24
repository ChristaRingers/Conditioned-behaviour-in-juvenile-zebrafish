%% Unit test for assesing whether time bins are correctly calculated
% Test function for the get_time_bins function in Anaylis_code_cr

% Main function:
function tests = analysis_code_test
tests = functiontests(localfunctions); % localfunctions is an automatic handle to all functions in this file
end

%% Local functions

% Setup function to pass variables needed for testing
function setupOnce(testCase)  % Set up one path for loading variables
% open a figure, for example
testCase.TestData.fn=FunObj;
testCase.TestData.origPath = 'test_variables.mat';
end

%Function 1 description:
%[xNew,yNew,dNew,heatmaps,vNew] = get_time_bins(X_vector,Y_vector,D_vector,V_vector,time_window,info_vector,min_size_bin);
function test_get_time_bins(testCase) % testcase is an object that is automatically generated for the purpose of testing

% Load all variables necessary for testing
load(testCase.TestData.origPath)

% Test for a random time_window
time_window = randi(10); % [in minutes] --> needs to be a multiple of the session time
[xNew,~,~,~,~] = testCase.TestData.fn.get_time_bins(X_vector,Y_vector,D_vector,V_vector,time_window,info_vector,min_size_bin); %#ok<USENS>

calculated_bins = size(xNew,1);

% Calculate an expected number of timebins
temp=0; expected_bins=0;
for t=1:size(D_vector,1)
    temp  = ceil(round(sum(D_vector{t,1})/(time_window*60000),1));
    expected_bins=expected_bins+temp;
end

verifyEqual(testCase,calculated_bins,expected_bins) % Same as assert in python
end


%Function 2 description:
%[distance] = calculate_distance(xpos, ypos);
function test_calculate_distance_horizontal(testCase)
start=1;step=randi(10);endpoint=1+step+step*randi(100);
expected_distance_horizontal=endpoint-start;
testHorizontal=start:step:endpoint;
ytest=ones(1,length(testHorizontal));

dist=testCase.TestData.fn.calculate_distance(testHorizontal,ytest);
verifyEqual(testCase,dist(1),step)
verifyEqual(testCase,sum(dist),expected_distance_horizontal)
end

function test_calculate_distance_vertical(testCase)
start=1;step=randi(10);endpoint=1+step+step*randi(100);
expected_distance_vertical=endpoint-start;
testVertical=start:step:endpoint;
xtest=ones(1,length(testVertical));

dist=testCase.TestData.fn.calculate_distance(xtest,testVertical);
verifyEqual(testCase,dist(1),step)
verifyEqual(testCase,sum(dist),expected_distance_vertical)
end


function test_calculate_distance_diagonal(testCase)
% xtest=randi([-10,10],1);
xtest=0.5-rand();
ytest=0.5-rand();

% ytest=randi([-10,10],1);
dist=testCase.TestData.fn.calculate_distance([0,xtest],[0,ytest]);
expected_value=abs(xtest/cos(atan(ytest/xtest)));
verifyEqual(testCase,expected_value,dist)
end

%% Test function template
% function testRealSolution(testCase)
% actSolution = quadraticSolver(1,-3,2);
% expSolution = [2 1];
% verifyEqual(testCase,actSolution,expSolution)
% end