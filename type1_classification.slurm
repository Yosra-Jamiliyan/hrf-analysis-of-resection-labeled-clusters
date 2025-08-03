#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --partition=cpu2019
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load fsl/6.0.4
export FSLDIR=~/fsl
export PATH=$FSLDIR/bin:$PATH
source $FSLDIR/etc/fslconf/fsl.sh

sub=40
typ=1

cluster_mask="/home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_cope_Ftest_sub${sub}.gfeat/cope1.feat/cluster_mask_zfstat${typ}.nii.gz"
zfstat="/home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_cope_Ftest_sub${sub}.gfeat/cope${typ}.feat/stats/zfstat${typ}.nii.gz"
pe1="/home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_cope_Ftest_sub${sub}.gfeat/cope${typ}.feat/stats/pe1.nii.gz"
pe2="/home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_cope_Ftest_sub${sub}.gfeat/cope${typ}.feat/stats/pe2.nii.gz"
pe3="/home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_cope_Ftest_sub${sub}.gfeat/cope${typ}.feat/stats/pe3.nii.gz"
resection_mask_orig="/home/yosra.jamiliyan1/new_data/ICE0${sub}/4_MRI/1_Anatomicals/PostOp/postop_mask_in_MNI.nii.gz"


cd /home/yosra.jamiliyan1/new_data/ICE0${sub}/custom_basis_fun/2nd_del_dis/

# Step 1 Resample resection mask to match zfstat.nii.gz file
#resampled_mask="resection_mask_resampled.nii.gz"
#flirt -in "$resection_mask_orig" -ref "$zfstat" -applyxfm -usesqform -out "$resampled_mask"

# Step 2 Make 2cm (10 voxel) expansion of resection mask
fslmaths "$resection_mask_orig" -kernel sphere 20 -dilF resection_mask_spherical_2cm.nii.gz

n_clusters=$(fslstats "$cluster_mask" -R | awk '{print int($2)}')

for i in $(seq 1 $n_clusters); do
  fslmaths "$cluster_mask" -thr $i -uthr $i -bin cluster_${i}.nii.gz
  fslmaths "$zfstat" -mas cluster_${i}.nii.gz z_in_cluster_${i}.nii.gz

  max_info=$(fslstats z_in_cluster_${i}.nii.gz -R -x)
  read -r _ _ X Y Z <<< "$max_info"
  X=$(printf "%.0f" "$X"); Y=$(printf "%.0f" "$Y"); Z=$(printf "%.0f" "$Z")

  fslmaths "$zfstat" -mul 0 -add 1 -roi $X 1 $Y 1 $Z 1 0 1 -bin max_voxel_mask_${i}.nii.gz

  vox_in_mask=$(fslstats "$resection_mask_orig" -k max_voxel_mask_${i}.nii.gz -V | awk '{print $1}')
  if (( vox_in_mask > 0 )); then
    label="Fully_Concordant"
  else
    vox_in_dilated=$(fslstats resection_mask_spherical_2cm.nii.gz -k max_voxel_mask_${i}.nii.gz -V | awk '{print $1}')
    if (( vox_in_dilated > 0 )); then
      label="Partially_Concordant"
    else
      overlap_cluster=$(fslstats "$resection_mask_orig" -k cluster_${i}.nii.gz -V | awk '{print $1}')
      if (( overlap_cluster > 0 )); then
        label="Partially_Discordant"
      else
        label="Fully_Discordant"
      fi
    fi
  fi

  echo "Cluster $i is $label"
  echo "Cluster $i is $label" > cluster_${i}_label.txt



  # Mask the PE image by the current cluster
  fslmaths "$pe1" -mas cluster_${i}.nii.gz pe1_cluster${i}.nii.gz
  fslmaths "$pe2" -mas cluster_${i}.nii.gz pe2_cluster${i}.nii.gz
  fslmaths "$pe3" -mas cluster_${i}.nii.gz pe3_cluster${i}.nii.gz


done
