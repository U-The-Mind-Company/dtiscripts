%% find_peak.m from Frank Yeh
%  http://dsi-studio.labsolver.org/Manual/Reconstruction
%  adapted for speed by Steven Baete, NYU SOM CBI, 2013-2017

% assumes odf_faces has been calculated as     
%     odf_faces = odf_faces + 1;
%     odf_faces = odf_faces - (odf_faces > length(odf))*length(odf);

    function p = find_peak_sb(odf,odf_faces)
    is_peak = odf;
    temp1 = odf(odf_faces(1,:));
    temp2 = odf(odf_faces(2,:));
    temp3 = odf(odf_faces(3,:));
    temp = temp2 >= temp1 | temp3 >= temp1;
    is_peak(odf_faces(1,temp)) = 0;
    temp = temp1 >= temp2 | temp3 >= temp2;
    is_peak(odf_faces(2,temp)) = 0;
    temp = temp2 >= temp3 | temp1 >= temp3;
    is_peak(odf_faces(3,temp)) = 0;
    %is_peak(odf_faces(1,odf(odf_faces(2,:)) >= odf(odf_faces(1,:)) | ...
    %    odf(odf_faces(3,:)) >= odf(odf_faces(1,:)))) = 0;
    %is_peak(odf_faces(2,odf(odf_faces(1,:)) >= odf(odf_faces(2,:)) | ...
    %    odf(odf_faces(3,:)) >= odf(odf_faces(2,:)))) = 0;
    %is_peak(odf_faces(3,odf(odf_faces(2,:)) >= odf(odf_faces(3,:)) | ...
    %    odf(odf_faces(1,:)) >= odf(odf_faces(3,:)))) = 0;
    [values,ordering] = sort(-is_peak);
    p = ordering(values < 0);
    end
