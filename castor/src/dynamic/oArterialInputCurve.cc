/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2021 all CASToR contributors listed below:

    --> Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Thibaut MERLIN, Mael MILLARDET, Simon STUTE, Valentin VIELZEUF, Zacharias CHALAMPALAKIS

This is CASToR version 3.1.1.
*/

/*!
  \file
  \ingroup  dynamic
  \brief    Implementation of class oArterialInputCurve
*/

#include "oArterialInputCurve.hh"
#include "sOutputManager.hh"
#include "gOptions.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn oArterialInputCurve
  \brief Constructor of oArterialInputCurve. Simply set all data members to default values.
*/

oArterialInputCurve::oArterialInputCurve()
{

  m_pathToAICfile = "";
  m_verbose = -1;

  m_nbinputDataPoints = -1;
  mp_AICDataPoints = NULL;
  mp_AICDataPointsTimes = NULL;


  mp_AICIntrpY = NULL;
  mp_AICIntrpY_downsampled = NULL;

  m_nbTimeFrames = -1;
  mp_frameTimeStartInMs = NULL;
  mp_frameTimeStopInMs = NULL;

}

oArterialInputCurve::~oArterialInputCurve()
{

  if (mp_AICDataPoints) delete(mp_AICDataPoints);
  if (mp_AICDataPointsTimes) delete(mp_AICDataPointsTimes);

  if (mp_AICIntrpY) delete(mp_AICIntrpY);

}

int oArterialInputCurve::CheckParameters()
{
  // Check that the provided sampled AIC covers the first frame
  if (mp_AICDataPointsTimes[0]>(mp_frameTimeStartInMs[0]))
  {
    Cerr("***** oArterialInputCurve::CheckParameters() -> AIC datapoints need to start before the first frame " << endl);
    return 1;
  }
  // Check if the provided sampled AIC covers the requested total reconstruction time ( last time point of last frame )
  if ( mp_frameTimeStopInMs[m_nbTimeFrames-1] > mp_AICDataPointsTimes[m_nbinputDataPoints-1] )
  {
    Cerr("***** oArterialInputCurve::CheckParameters() -> Requested reconstruction framing is longer than the provided AIC time points" << endl);
    return 1;
  }
  // End
  return 0;
}

void oArterialInputCurve::SetFrames(int a_nbTimeFrames, uint32_t *a_frameTimeStartInMs, uint32_t *a_frameTimeStopInMs)
{
  m_nbTimeFrames = a_nbTimeFrames;
  mp_frameTimeStartInMs = a_frameTimeStartInMs;
  mp_frameTimeStopInMs = a_frameTimeStopInMs;
}

int oArterialInputCurve::InitializeInputData()
{

  if (m_verbose>=2) Cout("oArterialInputCurve::InitializeInputData() -> Initializing input from the provided Arterial Input Curve " << endl );

  // Read number of AIC data points
  if (ReadDataASCIIFile(m_pathToAICfile,"AIC_number_of_points",&m_nbinputDataPoints,1,KEYWORD_MANDATORY))
  {
    Cerr("***** oArterialInputCurve::InitializeInputData() -> Error while reading number of data points from AIC file )" << endl);
    return 1;
  }

  // Allocate memory for data points
  mp_AICDataPoints = new HPFLTNB[m_nbinputDataPoints];
  mp_AICDataPointsTimes = new uint32_t [m_nbinputDataPoints];
  // need floating point numbers for input data points
  FLTNB* input_time_points_FLTNB = new FLTNB[ m_nbinputDataPoints ];


  // Populate arrays with data from AIC file
  if ( ReadDataASCIIFile(m_pathToAICfile, "AIC_data_points", mp_AICDataPoints, m_nbinputDataPoints , KEYWORD_MANDATORY) )
  {
    Cerr( "***** oArterialInputCurve::InitializeInputData() -> Error while trying to read Data points from AIC file: " << m_pathToAICfile <<  endl);
    return 1;
  }
  if ( ReadDataASCIIFile(m_pathToAICfile, "AIC_time_points", input_time_points_FLTNB, m_nbinputDataPoints , KEYWORD_MANDATORY) )
  {
    Cerr("***** oArterialInputCurve::InitializeInputData() -> Error while trying to read Time Data points from AIC file: " << m_pathToAICfile << endl);
    return 1;
  }

  // Get Units from file
  string Units = "";
  ReadDataASCIIFile(m_pathToAICfile, "AIC_units", &Units,1, KEYWORD_OPTIONAL);

  bool time_points_inMinutes = false;

  if (Units=="minutes" || Units=="Minutes")
  {
    time_points_inMinutes = true ;
  }

  // Case of Time points provided in minutes
  if (time_points_inMinutes)
  {
    // Convert all time points from minutes to milliseconds and cast into the uint32_t variable
    for (int i=0;i<m_nbinputDataPoints;i++) { mp_AICDataPointsTimes[i] = (uint32_t) (input_time_points_FLTNB[i]*60000);}
  }
  // Else case of time points provided in seconds
  else
  {
    // Convert all time points from seconds to milliseconds and cast into the uint32_t variable
    for (int i = 0; i < m_nbinputDataPoints; i++) { mp_AICDataPointsTimes[i] = (uint32_t) (input_time_points_FLTNB[i]*1000);}
  }

  // delete temporary array
  if (input_time_points_FLTNB) delete [] (input_time_points_FLTNB);

  // Print data points
  if (m_verbose>=3)
  {
    Cout("oArterialInputCurve::InitializeInputData() -> Printing input Data \n Time(ms),Datapoint (Bq/ml) " << endl );
    for (int i = 0; i < m_nbinputDataPoints; i++)
    {
        Cout (mp_AICDataPointsTimes[i] << ", " << mp_AICDataPoints [i]<< endl);
    }
  }
  // End
  return 0;
}

int oArterialInputCurve::InterpolateAIC()
{
  // Creating variables for storing the slopes and intercepts between each pair of datapoints
  HPFLTNB* slopes =  new HPFLTNB[m_nbinputDataPoints-1];
  HPFLTNB* intercepts = new HPFLTNB[m_nbinputDataPoints-1];

  // Allocating memory for the interpolated data values ( one for each discrete unit -> 1 millisecond )
  // Only calculating and storing the interpolated values until the end of the last frame
  mp_AICIntrpY = new HPFLTNB[mp_frameTimeStopInMs[m_nbTimeFrames-1]];

  // Calculating slopes and intercepts for each pair  & interpolate for every millisecond value
  for (int i=0;i<(m_nbinputDataPoints-1);i++)
  {
    // m = (y2-y1)/(x2-x1)
    slopes[i] = (mp_AICDataPoints[i + 1] - mp_AICDataPoints[i]) / ((HPFLTNB) (mp_AICDataPointsTimes[i + 1] - mp_AICDataPointsTimes[i]));
    //b = (-1)*(m)*(x2)+y2
    intercepts[i] = -1 * slopes[i] * ((HPFLTNB)(mp_AICDataPointsTimes[i + 1])) + mp_AICDataPoints[i + 1];

    // Calculate the interpolated values for the evaluated data-pair and populate mp_AICIntrpY
    // If first data point is within the range proceed
    if (mp_AICDataPointsTimes[i] <= mp_frameTimeStopInMs[m_nbTimeFrames - 1])
    {
      // Check if second data point is also within the range
      if (mp_AICDataPointsTimes[i + 1] <= mp_frameTimeStopInMs[m_nbTimeFrames - 1])
      {
        // Populate the whole range with interpolated values
        for (uint32_t k = mp_AICDataPointsTimes[i]; k < mp_AICDataPointsTimes[i + 1]; k++)
        {
          mp_AICIntrpY[k] = (HPFLTNB)k * slopes[i] + intercepts[i];
        }
      }
      // if second data point is out of range , populate up until the last frame point
      else if (mp_AICDataPointsTimes[i + 1] > mp_frameTimeStopInMs[m_nbTimeFrames - 1])
      {
        // Populate the whole range with interpolated values
        for (uint32_t k = mp_AICDataPointsTimes[i]; k <= mp_frameTimeStopInMs[m_nbTimeFrames - 1]; k++)
        {
          mp_AICIntrpY[k] = (HPFLTNB)k * slopes[i] + intercepts[i];
        }
      }
    }
  }
  // Report any negative values and set to zero
  for (uint32_t k = 0; k <= mp_frameTimeStopInMs[m_nbTimeFrames - 1]; k++)
  {
    if ( mp_AICIntrpY[k] <0 )
    {
      mp_AICIntrpY[k] = 0 ;
      Cout (" Negative interpolated AIF value detected at: " << k <<"ms --> setting AIF to zero" << endl);
    }
  }

  // clear memory
  if (slopes) delete[] (slopes);
  if (intercepts) delete[] (intercepts);
  // Also clear input data points from memory , as we dont re-use them for interpolation
  // (that could change in the future with new models that might require other operations on the input datapoints)
  if (mp_AICDataPoints) delete [] (mp_AICDataPoints);
  if (mp_AICDataPointsTimes) delete [] (mp_AICDataPointsTimes);


  //Normal End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//  Downsampling of interpolated AIC for speeding-up convolution operations
// =====================================================================
int oArterialInputCurve::Downsample()
{
  // Number of total samples required
  int total_samples_forDownsample = (int)(((float)mp_frameTimeStopInMs[m_nbTimeFrames - 1] / 1000.) * 10);
  // plus one for index 0
  total_samples_forDownsample+=1;
  if (m_verbose>=3) Cout(" Downsampling to 100 msec time intervals " << endl);
  // Allocate memory for downsampled AIC
  mp_AICIntrpY_downsampled = new HPFLTNB[total_samples_forDownsample];
  // Seting the index and the zero value (first data-point)
  int index=0;
  mp_AICIntrpY_downsampled[0] = mp_AICIntrpY[0] ;
  for (int i=1;i<=int(mp_frameTimeStopInMs[m_nbTimeFrames - 1]);i++)
  {
    // If exact division by 100 , datapoint is a 0.1 of a second -> set to downsampled datapoints
    if (i%100==0)
    {
      index ++;
      mp_AICIntrpY_downsampled[index] = mp_AICIntrpY[i];
    }
  }
  if (m_verbose>=3) Cout( " Downsampling complete to total # of points : " << index << endl);

  // Normal End
  return 0;
}