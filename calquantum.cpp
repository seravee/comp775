#include "calquantum.h"



//------------------------------------------------------------------------------------------------------------
        //code added by Xiao Yang for computing the uniformity
        //
        //
        //
        //what we need is 14 eigenvalues and from that form the 14 sigmas. the 14 eigenvalues
        //forms from 7 pairs:
        //1. area of medial quad and upper surface;
        //2. area of medial quad and lower surface;
        //3. length of vertical length of medial and upper;
        //4. length of vertical length of medial and lower;
        //5. length of horizental length of medial and upper;
        //6. length of horizental length of medial and lower;
        //7. distance between points on the edge of the medial sheet and on the crest

int calquantum(vtkSmartPointer<vtkSRep> srepfig, vtkSmartPointer<vtkSRepInterpolateMedialSheet> medialsheetinterpolator, vtkSmartPointer<vtkSRep> srepcrest, vtkSmartPointer<vtkSRepInterpolateMedialSpokesHermite> medialspokesinterpolator)
{
     unsigned CellNum = srepfig->GetNumberOfCells();
        int interlevel = medialsheetinterpolator->GetInterpolationLevel();
        double areas[CellNum];
        double upperareas[CellNum];
        double lowerareas[CellNum];

        vector<double> medialhori, medialver, upperhori, upperver, lowerhori, lowerver;


        for(unsigned i = 0; i < CellNum; i++)
        {
                areas[i] = 0;
                upperareas[i] = 0;
                lowerareas[i] = 0;
        }

        double step = pow((double)2, (double)interlevel);
        for(unsigned i = 0; i < CellNum; i++)
        {
                //calculate the area of each cell
                for(double u = 0; u <= step-1; u++)
                {
                        for(double v = 0; v <= step-1; v++)
                        {
                                vnl_vector<double> point[4];
                                point[0] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, u/step, v/step);
                                point[1] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, (u+1)/step, v/step);
                                point[2] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, (u+1)/step, (v+1)/step);
                                point[3] = medialsheetinterpolator->GetInterpolatedPoint((vtkIdType)i, u/step, (v+1)/step);

                                double midline = sqrt(pow(point[0][0]-point[2][0], 2) + pow(point[0][1]-point[2][1], 2) + pow(point[0][2]-point[2][2], 2));
                                double line11= sqrt(pow(point[0][0]-point[1][0], 2) + pow(point[0][1]-point[1][1], 2) + pow(point[0][2]-point[1][2], 2));
                                medialhori.push_back(line11);
                                double line12= sqrt(pow(point[0][0]-point[3][0], 2) + pow(point[0][1]-point[3][1], 2) + pow(point[0][2]-point[3][2], 2));
                                medialver.push_back(line12);
                                double line21= sqrt(pow(point[2][0]-point[1][0], 2) + pow(point[2][1]-point[1][1], 2) + pow(point[2][2]-point[1][2], 2));
                                medialver.push_back(line21);
                                double line22= sqrt(pow(point[2][0]-point[3][0], 2) + pow(point[2][1]-point[3][1], 2) + pow(point[2][2]-point[3][2], 2));
                                medialhori.push_back(line22);
                                //compute size using Heron's formula
                                double p1 = (midline+line11+line12)/2;
                                double p2 = (midline+line21+line22)/2;

                                double size1 = sqrt(abs(p1*(p1-midline)*(p1-line11)*(p1-line12)));
                                double size2 = sqrt(abs(p2*(p2-midline)*(p2-line21)*(p2-line22)));
                                areas[i] += size1+size2;

                                //compute upper side of the surface
                                vnl_vector<double> upperpoint[4];
                                upperpoint[0] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, u/step, v/step);
                                upperpoint[1] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, (u+1)/step, v/step);
                                upperpoint[2] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, (u+1)/step, (v+1)/step);
                                upperpoint[3] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 0, u/step, (v+1)/step);
                                for(int pointnum = 0; pointnum < 4; pointnum++)
                                {
                                        for(int pointdim = 0; pointdim < 3; pointdim++)
                                        {
                                                upperpoint[pointnum][pointdim] += point[pointnum][pointdim];
                                        }
                                }
                                midline = sqrt(pow(upperpoint[0][0]-upperpoint[2][0], 2) + pow(upperpoint[0][1]-upperpoint[2][1], 2) + pow(upperpoint[0][2]-upperpoint[2][2], 2));
                                line11= sqrt(pow(upperpoint[0][0]-upperpoint[1][0], 2) + pow(upperpoint[0][1]-upperpoint[1][1], 2) + pow(upperpoint[0][2]-upperpoint[1][2], 2));
                                upperhori.push_back(line11);
                                line12= sqrt(pow(upperpoint[0][0]-upperpoint[3][0], 2) + pow(upperpoint[0][1]-upperpoint[3][1], 2) + pow(upperpoint[0][2]-upperpoint[3][2], 2));
                                upperver.push_back(line12);
                                line21= sqrt(pow(upperpoint[2][0]-upperpoint[1][0], 2) + pow(upperpoint[2][1]-upperpoint[1][1], 2) + pow(upperpoint[2][2]-upperpoint[1][2], 2));
                                upperhori.push_back(line21);
                                line22= sqrt(pow(upperpoint[2][0]-upperpoint[3][0], 2) + pow(upperpoint[2][1]-upperpoint[3][1], 2) + pow(upperpoint[2][2]-upperpoint[3][2], 2));
                                upperver.push_back(line22);

                                p1 = (midline+line11+line12)/2;
                                p2 = (midline+line21+line22)/2;

                                size1 = sqrt(abs(p1*(p1-midline)*(p1-line11)*(p1-line12)));
                                size2 = sqrt(abs(p2*(p2-midline)*(p2-line21)*(p2-line22)));
                                upperareas[i] = size1+size2;

                                //compute lower side of the surface
                                vnl_vector<double> lowerpoint[4];
                                lowerpoint[0] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, u/step, v/step);
                                lowerpoint[1] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, (u+1)/step, v/step);
                                lowerpoint[2] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, (u+1)/step, (v+1)/step);
                                lowerpoint[3] = medialspokesinterpolator->GetInterpolatedSpoke((vtkIdType)i, 1, u/step, (v+1)/step);
                                for(int pointnum = 0; pointnum < 4; pointnum++)
                                {
                                        for(int pointdim = 0; pointdim < 3; pointdim++)
                                        {
                                                lowerpoint[pointnum][pointdim] += point[pointnum][pointdim];
                                        }
                                }
                                midline = sqrt(pow(lowerpoint[0][0]-lowerpoint[2][0], 2) + pow(lowerpoint[0][1]-lowerpoint[2][1], 2) + pow(lowerpoint[0][2]-lowerpoint[2][2], 2));
                                line11= sqrt(pow(lowerpoint[0][0]-lowerpoint[1][0], 2) + pow(lowerpoint[0][1]-lowerpoint[1][1], 2) + pow(lowerpoint[0][2]-lowerpoint[1][2], 2));
                                lowerhori.push_back(line11);
                                line12= sqrt(pow(lowerpoint[0][0]-lowerpoint[3][0], 2) + pow(lowerpoint[0][1]-lowerpoint[3][1], 2) + pow(lowerpoint[0][2]-lowerpoint[3][2], 2));
                                lowerver.push_back(line12);
                                line21= sqrt(pow(lowerpoint[2][0]-lowerpoint[1][0], 2) + pow(lowerpoint[2][1]-lowerpoint[1][1], 2) + pow(lowerpoint[2][2]-lowerpoint[1][2], 2));
                                lowerhori.push_back(line21);
                                line22= sqrt(pow(lowerpoint[2][0]-lowerpoint[3][0], 2) + pow(lowerpoint[2][1]-lowerpoint[3][1], 2) + pow(lowerpoint[2][2]-lowerpoint[3][2], 2));
                                lowerhori.push_back(line22);

                                p1 = (midline+line11+line12)/2;
                                p2 = (midline+line21+line22)/2;

                                size1 = sqrt(abs(p1*(p1-midline)*(p1-line11)*(p1-line12)));
                                size2 = sqrt(abs(p2*(p2-midline)*(p2-line21)*(p2-line22)));
                                lowerareas[i] = size1+size2;

                        }
                }
        }




        double sigma[14];
        int num = 0;
        //substract the mean from the data
        double medialsheetmean = 0;
        double uppermean = 0;
        double lowermean = 0;
        for(int i = 0; i < CellNum; i++)
        {
                medialsheetmean += areas[i]/CellNum;
                uppermean += upperareas[i]/CellNum;
                lowermean += lowerareas[i]/CellNum;
        }
        for(int i = 0; i < CellNum; i++)
        {
                areas[i] -=medialsheetmean;
                upperareas[i] -= uppermean;
                lowerareas[i] -= lowermean;
        }

        //form the matrix
        double matrix[2][2];
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        //calculate feature 1
        for(int i = 0; i < CellNum; i++)
        {
                matrix[0][0] += areas[i] * areas[i] / (CellNum-1);
                matrix[0][1] += areas[i] * upperareas[i] / (CellNum-1);;
                matrix[1][0] += areas[i] * upperareas[i] / (CellNum-1);;
                matrix[1][1] += upperareas[i] * upperareas[i] / (CellNum-1);;
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[0] << ' ' << sigma[1] << endl;

        //calculate feature 2
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < CellNum; i++)
        {
                matrix[0][0] += areas[i] * areas[i] / (CellNum-1);
                matrix[0][1] += areas[i] * lowerareas[i] / (CellNum-1);;
                matrix[1][0] += areas[i] * lowerareas[i] / (CellNum-1);;
                matrix[1][1] += lowerareas[i] * lowerareas[i] / (CellNum-1);;
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[2] << ' ' << sigma[3] << endl;

        //preprocessing the lines for medial, upper and lower sheet for matrix calculation
        int linesize = medialver.size();
        int medialvermean = 0, medialhorimean = 0, uppervermean = 0, upperhorimean = 0, lowervermean = 0, lowerhorimean = 0;
        for(int i = 0; i < linesize; i++)
        {
                medialvermean += medialver[i];
                medialhorimean += medialhori[i];
                uppervermean += upperver[i];
                upperhorimean += upperhori[i];
                lowervermean += lowerver[i];
                lowerhorimean += lowerhori[i];
        }

        medialvermean /= linesize;
        medialhorimean /= linesize;
        uppervermean /= linesize;
        upperhorimean /= linesize;
        lowervermean /= linesize;
        lowerhorimean /= linesize;

        for(int i = 0; i < linesize; i++)
        {
                medialver[i] -= medialvermean;
                medialhori[i] -= medialhorimean;
                upperver[i] -= uppervermean;
                upperhori[i] -= upperhorimean;
                lowerver[i] -= lowervermean;
                lowerhori[i] -= lowerhorimean;
        }

        //calculate feature 3
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < linesize; i++)
        {
                matrix[0][0] += medialver[i] * medialver[i] / (linesize - 1);
                matrix[0][1] += medialver[i] * upperver[i] / (linesize - 1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += upperver[i] * upperver[i] / (linesize - 1);
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[4] << ' ' << sigma[5] << endl;

        //calculate feature 4
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < linesize; i++)
        {
                matrix[0][0] += medialver[i] * medialver[i] / (linesize - 1);
                matrix[0][1] += medialver[i] * lowerver[i] / (linesize - 1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += lowerver[i] * lowerver[i] / (linesize - 1);
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[6] << ' ' << sigma[7] << endl;

        //calculate feature 5
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < linesize; i++)
        {
                matrix[0][0] += medialhori[i] * medialhori[i] / (linesize - 1);
                matrix[0][1] += medialhori[i] * upperhori[i] / (linesize - 1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += upperhori[i] * upperhori[i] / (linesize - 1);
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[8] << ' ' << sigma[9] << endl;


        //calculate feature 6
        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;
        for(int i = 0; i < linesize; i++)
        {
                matrix[0][0] += medialhori[i] * medialhori[i] / (linesize - 1);
                matrix[0][1] += medialhori[i] * lowerhori[i] / (linesize - 1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += lowerhori[i] * lowerhori[i] / (linesize - 1);
        }

        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[10] << ' ' << sigma[11] << endl;



        //calculate feature 7

        unsigned crestpointnum = srepcrest->GetNumberOfPoints();

        double pointdistances[crestpointnum];
        double crestpointdistances[crestpointnum];

        //compute the deviation of distance between interpolated points on the spoke
        for(unsigned i = 0; i < crestpointnum; i++)
        {
                unsigned next = (i+1)%crestpointnum;
                double point[3], nextpoint[3];
                srepcrest->GetPoint(i, point);
                srepcrest->GetPoint(next, nextpoint);
                pointdistances[i] = sqrt(pow(nextpoint[0]-point[0], 2) + pow(nextpoint[1]-point[1], 2) + pow(nextpoint[2] - point[2], 2));

                vtkSRep::VectorVNLType currentspokes = srepcrest->GetSpokes(i);
                vtkSRep::VectorVNLType nextspokes = srepcrest->GetSpokes(next);
                vtkSRep::VectorDoubleType radius = srepcrest->GetSpokesRadius(i);
                unsigned spokenum = currentspokes.size()/2;
                vtkSRep::VNLType p = currentspokes[spokenum] * radius[spokenum];
                for(int num = 0; num < 3; num++)
                {
                        point[num] += p[num];
                }

                radius = srepcrest->GetSpokesRadius(next);
                p = nextspokes[spokenum] * radius[spokenum];
                for(int num = 0; num < 3; num++)
                {
                        nextpoint[num] += p[num];
                }
                crestpointdistances[i] = sqrt(pow(nextpoint[0]-point[0], 2) + pow(nextpoint[1]-point[1], 2) + pow(nextpoint[2] - point[2], 2));
        }

        double pointsmean = 0;
        double crestpointsmean= 0;

        for(unsigned i = 0; i < crestpointnum; i++)
        {
                pointsmean += pointdistances[i];
                crestpointsmean += crestpointdistances[i];
        }
        pointsmean /= crestpointnum;
        crestpointsmean /= crestpointnum;

        for(unsigned i = 0; i < crestpointnum; i++)
        {
                pointdistances[i] -= pointsmean;
                crestpointdistances[i] -= crestpointsmean;
        }

        matrix[0][0] = matrix[0][1] = matrix[1][0] = matrix[1][1] = 0;

        for(unsigned i = 0; i < crestpointnum; i++)
        {
                matrix[0][0] += pointdistances[i] * pointdistances[i] / (crestpointnum-1);
                matrix[0][1] += pointdistances[i] * crestpointdistances[i] / (crestpointnum-1);
                matrix[1][0] = matrix[0][1];
                matrix[1][1] += crestpointdistances[i] * crestpointdistances[i] / (crestpointnum-1);
        }
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) + sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        sigma[num++] = ((matrix[0][0] + matrix[1][1]) - sqrt( pow(matrix[0][0]-matrix[1][1], 2) + 4 * pow(matrix[0][1], 2)))/2;
        cout << sigma[12] << ' ' << sigma[13] << endl;

        double feature_entropy = 0;

        for(int i = 0; i < 14; i++)
        {
                feature_entropy += log(sigma[i]);
        }
        feature_entropy += 14 * (0.5 + 0.5 * log(3.1415926535));


        return feature_entropy;
}
