for file in $(echo sub-CON001_run-0_eyes_only_mc); #$(cat only_mc_list.txt);
do 
	file_full=$(echo $file | cut -d '.' -f 1)
	echo "$file_full"

	3dBlurToFWHM -input only_mc/"$file_full".nii.gz -FWHM 6 -prefix only_mc/smoothed/"$file_full"_sm.nii.gz -overwrite

	for side in right left;
	do 
		3dcalc -a ../nau_mask_"$side".nii.gz -b only_mc/smoothed/"$file_full"_sm.nii.gz -expr '(a*b)' -prefix only_mc/separate_eyes/"$file_full"_"$side".nii.gz -overwrite


		3dvolreg -input only_mc/separate_eyes/"$file_full"_"$side".nii.gz -base only_mc/separate_eyes/"$file_full"_"$side".nii.gz[0] -1Dfile only_mc/within_mask_mc/"$file_full"_mc_"$side".txt -nomaxdisp -prefix only_mc/within_mask_mc/"$file_full"_mc_"$side".nii.gz -overwrite
	done
done
