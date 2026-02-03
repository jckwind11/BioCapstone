main("SC");
main("FC");

function [fmtStr] = genFmtStr(nc)
    res = "";
    for i = 1:nc
        res = append(res, "%1.4f,");
    end
    fmtStr = extractBetween(res, 1, strlength(res)-1) + "\n";
end       

function [] = main(type)
    file = load('Individual_Connectomes.mat');
    fb = "Individual_Connectomes_";
    sizes = ["68", "114", "219", "448", "1000"];

    if (strcmp(type, "SC"))
        fn = append(fb, "SC_");
        data = file.connMatrices.SC(:,1);
    else 
        fn = append(fb, "FC_");
        data = file.connMatrices.FC(:,1);
    end

    for n = 1:length(sizes)
        C = data{n,:,:,:};
        [~, ncols, nsubjects] = size(C);
    
        fileID = fopen(fn+sizes(n)+'.csv','w');
        [fmtStr] = genFmtStr(ncols);
    
        for subject = 1:nsubjects
            fprintf(fileID, fmtStr, C(:,:,subject));
        end
        
        disp(type+" "+sizes(n)+" done.");
        fclose(fileID);
    end
end 

