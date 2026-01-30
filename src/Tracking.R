args<-commandArgs(TRUE)
ex <-args[1] 
print(ex)
ds<-as.character(args[2]) 
de<-as.character(args[3]) 
store.folder<-as.character(args[4]) 
print(paste(ex,ds,de,store.folder))


# -----------------------
#ex<-"SWOT_WATER"
#ds<- as.character("20230801")
#de<- as.character("20230803")
#store.folder<-"/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER/SWOT_WATER_LOW_CY_OUT"
#print(paste(ex,ds,de,store.folder))
# -----------------------


MYD<- sprintf("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER")
# Namelist
namelist<-sprintf("%s/NAMELIST_%s.R",MYD,ex)
print(namelist)
source(namelist)

#ex<-"SSH_DYN"
#ds<-"20230330"
#de<-"20230710"
#store.folder<-sprintf("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY/OUT")


.libPaths("/gpfs/home/bonaduce/R_library")
library(ncdf4,lib.loc="/gpfs/home/bonaduce/R_library")
library(SDMTools,lib.loc="/gpfs/home/bonaduce/R_library")
library(fields,lib.loc="/gpfs/home/bonaduce/R_library")
library(date,lib.loc="/gpfs/home/bonaduce/R_library")

starting_date<-ds # as.integer(args[2]) ## first frame read
ending_date<-de # as.integer(args[3])  ## last frame read

startdate<-as.Date(starting_date,format="%Y%m%d")
enddate<-as.Date(ending_date,format="%Y%m%d")
#number of total frame

# Check time-window
command<-sprintf("ls -1 %s/*_%s.nc",dirONE,ex_lab)
listfile<-system(command,inter=TRUE)
listfile<-basename(listfile)
npiece<-length(unlist(strsplit(listfile[1],split="_")))
list_days<-unlist(lapply(strsplit(listfile,split="_"), "[", 1))
list_days<-unlist(lapply(strsplit(list_days,split="\\."), "[", 1))
list_days<-as.Date(list_days,format="%Y%m%d")
idx_days <- list_days >= startdate & list_days <= enddate
list_days <- list_days[idx_days]
nd <- length(list_days)
startdate <- list_days[1]
enddate <- list_days[nd]
step_time <- as.numeric(list_days[2] - list_days[1])  

n_times<- nd

source(sprintf("%s/fuction_tt2.R",MYD))
source(sprintf("%s/Ini.R",MYD))
exp.folder <- sprintf("%s/",store.folder)
print(exp.folder)
if(!file.exists(exp.folder)) { dir.create(exp.folder) }

#max_number of eddies
max_eddies=1000
#max dist in 1day
DIST<-100*1000 #20 km
#max dist west in 1day
DIST_west<-100*1000 #20 km

#for boundary purpouse
max_dist<-1
if (ex == "ALT_GLO") {
 max_dist<-11
}

#for selecting eddies with similar shape
ncell_insteadof_dist_fraction_max<-10 #2times bigger
ncell_insteadof_dist_fraction_min<-0.1

ncell_fraction_max_cutout<-225# 5times bigger
ncell_fraction_min_cutout<-0.0020


#source(file="filled_contour2.R")

## files
file_median_x<-sprintf("%s/median_x.dat",exp.folder)
file_median_y<-sprintf("%s/median_y.dat",exp.folder)
file_ncell<-sprintf("%s/numcell.dat",exp.folder)
file_clockwise_eddy<-sprintf("%s/clockwise_eddy.dat",exp.folder)
file_eddy_index<-sprintf("%s/eddy_index.dat",exp.folder)

filemask<-sprintf("%s/meshmask.nc",dirTWO)

ncm<-nc_open(filemask)
#lon<-ncm$dim$lon$vals
#lat<-ncm$dim$latvals
depth_value<-1 #ncm$dim$depth$vals
#depth_value<-depth_value[1:z_window]

nav_lon<-ncvar_get(ncm,"lon")
nav_lat<-ncvar_get(ncm,"lat")
#dx<-ncvar_get(ncm,"e1t")
#dy<-ncvar_get(ncm,"e2t")
nc_close(ncm)
dz<-1

matrix_lon<-nav_lon
matrix_lat<-nav_lat

boardcond<-matrix_lon[dim(nav_lon)[1],]
matrix_lon_up<-my_shift.up(matrix_lon,1)
matrix_lon_up[1,]<-boardcond

boardcond<-rev(matrix_lat[,dim(nav_lat)[2]]-1)
matrix_lat_left<-my_shift.left(matrix_lat,1)
matrix_lat_left[,dim(nav_lat)[2]]<-boardcond

source(sprintf("%s/hav.R",MYD))
dy<-haversine(matrix_lon,matrix_lon,matrix_lat,matrix_lat_left)*1000
dx<-haversine(matrix_lon,matrix_lon_up,matrix_lat,matrix_lat)*1000
print(dim(dy))
dx[1,]<-dx[2,]
dx[dim(dx)[1],]<-dx[dim(dx)[1]-1,]
dy[,1]<-dx[,2]
print(dim(dy)[1]-1)
#dy[,dim(dy)[2]]<-dy[,dim(dy)[1]-1]
dy[,dim(dy)[2]]<-dy[,dim(dy)[2]-1]

nav_dx<-dx
nav_dy<-dy

rm(matrix_lon,matrix_lat,matrix_lon_up,matrix_lat_left,boardcond)

x_window<-dim(nav_lon)[1]
y_window<-dim(nav_lat)[2]
print(c(x_window,y_window))

# converting distance in number of grid points
max_DIST_x<-DIST/nav_dx
max_DIST_y<-DIST/nav_dy

# converting distance towards west in number of grid points
max_DIST_x_west<-DIST_west/nav_dx
max_DIST_y_west<-DIST_west/nav_dy

# tot num of eddies per frame
eddies_tot_num<-array(NA,dim=c(n_times))

#
# reading files
#
row_data_median_x<-scan(file_median_x)
row_data_median_y<-scan(file_median_y)
row_data_ncell<-scan(file_ncell)
row_data_clockwise_eddy<-scan(file_clockwise_eddy)

#--- algorithm for reading files
data_n<-0
for (i in 1:n_times)
{
data_n<-data_n+1
eddies_tot_num[i]<-row_data_median_x[data_n]
data_n<-data_n+eddies_tot_num[i]
}

max_eddies<-max(eddies_tot_num)
#max_eddies<-max(eddies_tot_num,na.rm=TRUE)
print(max_eddies)
#array storing input data
median_x_t<-array(NA,dim=c(n_times,max_eddies))
median_y_t<-array(NA,dim=c(n_times,max_eddies))
ncell_t<-array(NA,dim=c(n_times,max_eddies))
index_ed_t<-array(NA,dim=c(n_times,max_eddies))
cycl_eddy_t<-array(NA,dim=c(n_times,max_eddies))

data_n<-0
for (i in 1:n_times)
{
    data_n<-data_n+1
    if(eddies_tot_num[i]>0)
    {
    median_x_t[i,1:eddies_tot_num[i]]<-row_data_median_x[(data_n+1):(data_n+eddies_tot_num[i])]
    median_y_t[i,1:eddies_tot_num[i]]<-row_data_median_y[(data_n+1):(data_n+eddies_tot_num[i])]
    ncell_t[i,1:eddies_tot_num[i]]<-row_data_ncell[(data_n+1):(data_n+eddies_tot_num[i])]
    cycl_eddy_t[i,1:eddies_tot_num[i]]<-row_data_clockwise_eddy[(data_n+1):(data_n+eddies_tot_num[i])]
    }
    data_n<-data_n+eddies_tot_num[i]
}

rm(row_data_median_x,row_data_median_y,row_data_ncell,row_data_clockwise_eddy)

#
# files now read and stored in above files
#

#relative indexes
rel_time<-0
max_id_prev_edd<-0
for (rel_time in seq(1, n_times, by=1))
    {
    start.time<-Sys.time()

    #each eddy in any frame has its own number
    if(rel_time>1) max_id_prev_edd<-max_id_prev_edd+eddies_tot_num[rel_time-1]
    if(eddies_tot_num[rel_time]>0)  index_ed_t[rel_time,1:eddies_tot_num[rel_time]]<-1:eddies_tot_num[rel_time]+max_id_prev_edd
    #eddies of the first record are independent by definition
    if (rel_time==1)
        {
		write(eddies_tot_num[rel_time], file = file_eddy_index,ncolumns = 1, sep = " ")
        write(index_ed_t[rel_time,1:eddies_tot_num[rel_time]], file = file_eddy_index,ncolumns = eddies_tot_num[rel_time],append = TRUE, sep = " ")
        next
        }
    cat("rel time=",rel_time)
    if(eddies_tot_num[rel_time-1]>0)
    {
    prev_pos<- array(cbind(median_x_t[rel_time-1,1:eddies_tot_num[rel_time-1]],median_y_t[rel_time-1,1:eddies_tot_num[rel_time-1]]),dim=c(eddies_tot_num[rel_time-1],2))
    prev_ncell<-array(ncell_t[rel_time-1,1:eddies_tot_num[rel_time-1]],dim=eddies_tot_num[rel_time-1])
    prev_cycl<-cycl_eddy_t[rel_time-1,1:eddies_tot_num[rel_time-1]]
    }
    if(eddies_tot_num[rel_time]>0)
    {
    new_pos <-array(cbind(median_x_t[rel_time,1:eddies_tot_num[rel_time]],median_y_t[rel_time,1:eddies_tot_num[rel_time]]),dim=c(eddies_tot_num[rel_time],2))
    new_ncell <-array(ncell_t[rel_time,1:eddies_tot_num[rel_time]],dim=eddies_tot_num[rel_time])
    new_cycl<-cycl_eddy_t[rel_time,1:eddies_tot_num[rel_time]]
    }
    #boundary treatment
    if (eddies_tot_num[rel_time]>0 & eddies_tot_num[rel_time-1]>0 )
    {
    if(eddies_tot_num[rel_time-1]<eddies_tot_num[rel_time])
            {
            #-----MIRROR EDDIES :create mirror eddies and add om the right
            mirror_eddy_cond<- new_pos[,1]<=max_dist # find mirror eddies
            #--- LINK beetween mirror(2nd column) and mirrored(1st column)
            mirror_col<-which(mirror_eddy_cond)     # indexes on the matrix
            mirror_col<-cbind(mirror_col,(eddies_tot_num[rel_time]+1):(eddies_tot_num[rel_time]+length(mirror_col)))# related index on new_pos_with_mirror
    
            #---create extended positions , a.m., ncell arrays with mirror
            mirror_eddy_pos<-array(new_pos[mirror_eddy_cond,],dim=c((length(new_pos[mirror_eddy_cond,]) %/% 2),2)) # x,y position
            mirror_eddy_pos[,1]<- mirror_eddy_pos[,1]+x_window#only x_window we have already cut previous term
            mirror_new_cycl<-new_cycl[mirror_eddy_cond]
            mirror_new_ncell<-new_ncell[mirror_eddy_cond]
    
            #---
            new_pos_with_mirror<-rbind(new_pos,mirror_eddy_pos)
            new_cycl_withmirror<-c(new_cycl,mirror_new_cycl)
            new_ncell_with_mirror<-c(new_ncell,mirror_new_ncell)
    
    
            #----calculating dist matrix
            prev_matr_x<-array(prev_pos[,1],dim=c(length(prev_pos[,1]),length(new_pos_with_mirror[,1])))
            prev_matr_y<-array(prev_pos[,2],dim=c(length(prev_pos[,2]),length(new_pos_with_mirror[,2])))
            new_pos_with_mirror_matr_x<-t(array(new_pos_with_mirror[,1],dim=c(length(new_pos_with_mirror[,1]),length(prev_pos[,1]))))
            new_pos_with_mirror_matr_y<-t(array(new_pos_with_mirror[,2],dim=c(length(new_pos_with_mirror[,2]),length(prev_pos[,2]))))
    
            #instead of real distance we use the "infinite" distance with abs, differnt definition of distance are equivalent in R2
            dist_matr_x<-abs(prev_matr_x-new_pos_with_mirror_matr_x)
            dist_matr_y<-abs(prev_matr_y-new_pos_with_mirror_matr_y)
    
            round_prev_pos<-round(prev_pos)
            MAX_xdist_prev_west<-max_DIST_x_west[round_prev_pos]
            MAX_ydist_prev_west<-max_DIST_y_west[round_prev_pos]
            MAX_xdist_prev_east_north<-max_DIST_x[round_prev_pos]
            MAX_ydist_prev_east_north<-max_DIST_y[round_prev_pos]
            MAX_xdist_matrix_west<-array(MAX_xdist_prev_west,dim(dist_matr_x))
            MAX_ydist_matrix_west<-array(MAX_ydist_prev_west,dim(dist_matr_y))
    
            MAX_xdist_matrix_east_north<-array(MAX_xdist_prev_east_north,dim(dist_matr_x))
            MAX_ydist_matrix_east_north<-array(MAX_ydist_prev_east_north,dim(dist_matr_y))
    
            dist_cond<-dist_matr_x<MAX_xdist_matrix_west & dist_matr_y<MAX_ydist_matrix_east_north
            #---filtering
            #first filter on distance
            prev_new_within_dist<-which(dist_cond,arr.ind=TRUE)
            prev_new_within_dist<-array(prev_new_within_dist,dim=c(length(prev_new_within_dist)/2,2))
            
            prev_new_within_dist_x<-cbind(prev_pos[prev_new_within_dist[,1],1],new_pos_with_mirror[prev_new_within_dist[,2],1])
    
            cond_east_propapation<-(prev_new_within_dist_x[,2]-prev_new_within_dist_x[,1])<MAX_xdist_matrix_east_north[prev_new_within_dist] | (prev_new_within_dist_x[,2]-prev_new_within_dist_x[,1])>x_window/2
    
            prev_new_within_dist<-prev_new_within_dist[cond_east_propapation,]
            prev_new_within_dist<-array(prev_new_within_dist,dim=c(length(prev_new_within_dist)/2,2))
            #first filter on  a.m.
            prev_new_within_dist<-prev_new_within_dist[(abs(prev_cycl[prev_new_within_dist[,1]] - new_cycl_withmirror[prev_new_within_dist[,2]])<2),]
            prev_new_within_dist<-array(prev_new_within_dist,dim=c(length(prev_new_within_dist)/2,2))
            
            new_ncell_with_mirror_filtered<-new_ncell_with_mirror[prev_new_within_dist[,2]]
            prev_ncell_filtered<-prev_ncell[prev_new_within_dist[,1]]
    
            #second filter on number on cells and distance
            ncell_difference_filtered<- abs(new_ncell_with_mirror_filtered-prev_ncell_filtered)
            dist_filtered<-sqrt(dist_matr_x[prev_new_within_dist]^2+dist_matr_y[prev_new_within_dist]^2)#dist_matr[prev_new_within_dist]#
        if (length(ncell_difference_filtered)>0)
        {
            absol_indexes_filtered<-1:length(ncell_difference_filtered)
            idf<-factor(prev_new_within_dist[,1])
            level<-as.integer(levels(idf))
    
            # find index with minimum in dist and ncel difference
            rel_pos_min_dist<-array(tapply(dist_filtered,idf,which.min))
            rel_pos_min_ncell<-array(tapply(ncell_difference_filtered,idf,which.min))
    
            # listing for finding absol index and playing on rel pos
            ncell_difference_filtered_list<-array(tapply(ncell_difference_filtered,idf,c))
            prev_ncell_filtered_list<-array(tapply(prev_ncell_filtered,idf,c))# should be always the same value
            new_ncell_filtered_list<-array(tapply(new_ncell_with_mirror_filtered,idf,c))
    
            dist_filtered_list<-array(tapply(dist_filtered,idf,c))
            absol_indexes_filtered_list<-array(tapply(absol_indexes_filtered,idf,c))
    
            final_absol_index<-array(0,dim=c(length(level)))
    
            for( ab in 1:length(level))
                {
                temp_dist<-unlist(dist_filtered_list[ab])
                temp_diffncell<-unlist(ncell_difference_filtered_list[ab])
                temp_newncell<-unlist(new_ncell_filtered_list[ab])
                temp_prevncell<-unlist(prev_ncell_filtered_list[ab])#should be aways the same value
                temp_prevncell<-temp_prevncell[1]
                temp_perc_ncell<-temp_newncell/temp_prevncell
                #starting with position of min_dist
                rel_min<-rel_pos_min_dist[ab]
                # changing this minimum eventually
                if ( rel_pos_min_dist[ab] != rel_pos_min_ncell[ab])
                    {
                    dist_mincell<-temp_dist[rel_pos_min_ncell[ab]]
                    dist_mindist<-temp_dist[rel_pos_min_dist[ab]]
            
                    if ( abs(dist_mincell-dist_mindist)<6 & (temp_perc_ncell[rel_pos_min_dist[ab]]>ncell_insteadof_dist_fraction_max | temp_perc_ncell[rel_pos_min_dist[ab]]<ncell_insteadof_dist_fraction_min))    rel_min<-rel_pos_min_ncell[ab]
                        
                    if (abs(dist_mincell-dist_mindist)<1) rel_min<-rel_pos_min_ncell[ab]
                        
                    }
                #find absolute index of chosen rel_min
                temp_absol_pos<-unlist(absol_indexes_filtered_list[ab])
                final_absol_index[ab]<-temp_absol_pos[rel_min]
                if (temp_perc_ncell[rel_min]>ncell_fraction_max_cutout | temp_perc_ncell[rel_min]<ncell_fraction_min_cutout) final_absol_index[ab]<-NA
                
                }
                
            final_absol_index<-final_absol_index[!is.na(final_absol_index)]
            prev_new_within_dist<-prev_new_within_dist[final_absol_index,]
            prev_new_within_dist<-array(prev_new_within_dist,dim=c(length(prev_new_within_dist)/2,2))
            final_prev_new_withmirror<-prev_new_within_dist
    
            #eliminating duplicated columns
            if (sum(duplicated(final_prev_new_withmirror[,2]))>0)
                {
                cond<-final_prev_new_withmirror[,2] %in% final_prev_new_withmirror[duplicated(final_prev_new_withmirror[,2]),2]
                duplic_new_eddies<-final_prev_new_withmirror[cond,]
                idf<-factor(duplic_new_eddies[,2])
                level<-as.integer(levels(idf))
                absol_index<-1:length(duplic_new_eddies[,2])
                dist_duplic_new_eddies<-sqrt(dist_matr_x[duplic_new_eddies]^2+dist_matr_y[duplic_new_eddies]^2)
                rel_min_dist_list<-array(tapply(dist_duplic_new_eddies,idf,which.min))
                absol_index_list<-array(tapply(absol_index,idf,c))
                final_absol_index<-array(0,dim=c(length(level)))
        
                for( ab in 1:length(level))
                    {
                    temporal<-unlist(absol_index_list[ab])
                    final_absol_index[ab]<-temporal[rel_min_dist_list[ab]]
                    }
                eddies_todelete_index<-1:dim(duplic_new_eddies)[1] #row are not duplicated
                eddies_to_delete<-duplic_new_eddies[!(eddies_todelete_index %in% final_absol_index),]
                eddies_to_delete<-array(eddies_to_delete,dim=c(length(eddies_to_delete)%/%2,2))
                final_prev_new_withmirror<-final_prev_new_withmirror[!(final_prev_new_withmirror[,1] %in% eddies_to_delete[,1]),]
                }
            #---final_prev_new_withmirror has got non duplicated eddies, now we need to discern mirror
    
            dim(final_prev_new_withmirror)<-c(length(final_prev_new_withmirror) %/% 2 ,2)
            final_prev_new<-final_prev_new_withmirror
            found_col_mirror<-final_prev_new_withmirror[final_prev_new[,2]>eddies_tot_num[rel_time],]
            found_col_mirror<-array(found_col_mirror,c((length(found_col_mirror) %/% 2),2))
            pos_mirror_in<-which(final_prev_new_withmirror[,2]>eddies_tot_num[rel_time])
    
            if (length(pos_mirror_in)==0) {rel_pos_i_b<-NULL}else { rel_pos_i_b<- 1:length(pos_mirror_in)}
            for (i_b in rel_pos_i_b)
                {
                col_mirrored<-mirror_col[ mirror_col[,2] == found_col_mirror[i_b,2],1]
        
                if (sum(final_prev_new_withmirror[,2]==col_mirrored)>0 )
                    {
                    found_col_mirrored<-final_prev_new_withmirror[final_prev_new_withmirror[,2]>eddies_tot_num[rel_time],]
                    found_col_mirrored<-array(found_col_mirrored,c((length(found_col_mirrored) %/% 2),2))
                    pos_mirrored_in<-which(final_prev_new_withmirror[,2]==col_mirrored)
                    min_mirror<-sqrt(dist_matr_x[found_col_mirror[i_b,1],found_col_mirror[i_b,2]]^2+dist_matr_y[found_col_mirror[i_b,1],found_col_mirror[i_b,2]]^2)
                    min_mirrored<-sqrt(dist_matr_x[found_col_mirrored[1],found_col_mirrored[2]]^2+dist_matr_y[found_col_mirrored[1],found_col_mirrored[2]]^2)
                    
                    if (min_mirror<min_mirrored)
                        {
                        final_prev_new[pos_mirror_in[i_b],2]<-col_mirrored
                        final_prev_new[pos_mirrored_in,]<-c(NA,NA)
                        }
                    else
                        {
                        final_prev_new[pos_mirror_in[i_b],]<-c(NA,NA)
                        }
                    }
                else
                    {
                    final_prev_new[pos_mirror_in[i_b],2]<-col_mirrored
                    }
                }
            final_prev_new<-final_prev_new[!is.na(final_prev_new[,1]),]
            dim(final_prev_new)<-c(length(final_prev_new) %/% 2 ,2)
            
            #final results
            index_ed_t[rel_time,final_prev_new[,2]]<-index_ed_t[rel_time-1,final_prev_new[,1]]
        }
            }
        else #the other case
            {
            #-----MIRROR EDDIES :create mirror eddies and add om the right
            mirror_eddy_cond<- prev_pos[,1]<=max_dist # find mirror eddies
            #--- LINK beetween mirror(2nd column) and mirrored(1st column)
            mirror_col<-which(mirror_eddy_cond)     # indexes on the matrix
            mirror_col<-cbind(mirror_col,(eddies_tot_num[rel_time-1]+1):(eddies_tot_num[rel_time-1]+length(mirror_col)))# related index on new_pos_with_mirror

            #---create extended positions , a.m., ncell arrays with mirror
            mirror_eddy_pos<-array(prev_pos[mirror_eddy_cond,],dim=c((length(prev_pos[mirror_eddy_cond,]) %/% 2),2)) # x,y position
            mirror_eddy_pos[,1]<- mirror_eddy_pos[,1]+x_window#only x_window we have already cut previous term
            mirror_prev_cycl<-prev_cycl[mirror_eddy_cond]
            mirror_prev_ncell<-prev_ncell[mirror_eddy_cond]
        
            prev_pos_with_mirror<-rbind(prev_pos,mirror_eddy_pos)
            prev_cycl_withmirror<-c(prev_cycl,mirror_prev_cycl)
            prev_ncell_with_mirror<-c(prev_ncell,mirror_prev_ncell)
            #----calculating dist matrix
            new_matr_x<-array(new_pos[,1],dim=c(length(new_pos[,1]),length(prev_pos_with_mirror[,1])))
            new_matr_y<-array(new_pos[,2],dim=c(length(new_pos[,2]),length(prev_pos_with_mirror[,2])))
            prev_pos_with_mirror_matr_x<-t(array(prev_pos_with_mirror[,1],dim=c(length(prev_pos_with_mirror[,1]),length(new_pos[,1]))))
            prev_pos_with_mirror_matr_y<-t(array(prev_pos_with_mirror[,2],dim=c(length(prev_pos_with_mirror[,2]),length(new_pos[,2]))))
        
            dist_matr_x<-abs(new_matr_x-prev_pos_with_mirror_matr_x)
            dist_matr_y<-abs(new_matr_y-prev_pos_with_mirror_matr_y)
        
            round_new_pos<-round(new_pos)
            MAX_xdist_new_east_north<-max_DIST_x[round_new_pos]
            MAX_ydist_new_east_north<-max_DIST_y[round_new_pos]
            MAX_xdist_new_west<-max_DIST_x_west[round_new_pos]
            MAX_ydist_new_west<-max_DIST_y_west[round_new_pos]
        
            MAX_xdist_matrix_west<-array(MAX_xdist_new_west,dim(dist_matr_x))
            MAX_ydist_matrix_west<-array(MAX_ydist_new_west,dim(dist_matr_y))
            MAX_xdist_matrix_east_north<-array(MAX_xdist_new_east_north,dim(dist_matr_x))
            MAX_ydist_matrix_east_north<-array(MAX_ydist_new_east_north,dim(dist_matr_y))
            
            dist_cond<-dist_matr_x<MAX_xdist_matrix_west & dist_matr_y<MAX_ydist_matrix_east_north
            #---filtering
            #first filter on distance
            new_prev_within_dist<-which(dist_cond,arr.ind=TRUE)
            new_prev_within_dist<-array(new_prev_within_dist,dim=c(length(new_prev_within_dist)/2,2))
            new_prev_within_dist_x<-cbind(new_pos[new_prev_within_dist[,1],1], prev_pos_with_mirror[new_prev_within_dist[,2],1])
            
            cond_east_propapation<-(new_prev_within_dist_x[,1]-new_prev_within_dist_x[,2])<MAX_xdist_matrix_east_north[new_prev_within_dist] | (new_prev_within_dist_x[,1]-new_prev_within_dist_x[,2])>x_window/2
        
            new_prev_within_dist<-new_prev_within_dist[cond_east_propapation,]
            new_prev_within_dist<-array(new_prev_within_dist,dim=c(length(new_prev_within_dist)/2,2))
            #on a.m.
            new_prev_within_dist<-new_prev_within_dist[(abs(new_cycl[new_prev_within_dist[,1]] - prev_cycl_withmirror[new_prev_within_dist[,2]])<2),]
            new_prev_within_dist<-array(new_prev_within_dist,dim=c(length(new_prev_within_dist)/2,2))
            prev_ncell_with_mirror_filtered<-prev_ncell_with_mirror[new_prev_within_dist[,2]]
            new_ncell_filtered<-new_ncell[new_prev_within_dist[,1]]
            #second filter on number on cells and distance
            ncell_difference_filtered<- abs(prev_ncell_with_mirror_filtered-new_ncell_filtered)
            dist_filtered<-sqrt(dist_matr_x[new_prev_within_dist]^2+dist_matr_y[new_prev_within_dist]^2)
        if (length(ncell_difference_filtered)>0)
        {
            absol_indexes_filtered<-1:length(ncell_difference_filtered)
            idf<-factor(new_prev_within_dist[,1])
            level<-as.integer(levels(idf))
        
            # find index with minimum in dist and ncel difference
            rel_pos_min_dist<-array(tapply(dist_filtered,idf,which.min))
            rel_pos_min_ncell<-array(tapply(ncell_difference_filtered,idf,which.min))
        
            # listing for finding absol index and playing on rel pos
            ncell_difference_filtered_list<-array(tapply(ncell_difference_filtered,idf,c))
            new_ncell_filtered_list<-array(tapply(new_ncell_filtered,idf,c))# should be always the same value
            prev_ncell_filtered_list<-array(tapply(prev_ncell_with_mirror_filtered,idf,c))
        
            dist_filtered_list<-array(tapply(dist_filtered,idf,c))
            absol_indexes_filtered_list<-array(tapply(absol_indexes_filtered,idf,c))
            final_absol_index<-array(0,dim=c(length(level)))
        
            for( ab in 1:length(level))
                {
                temp_dist<-unlist(dist_filtered_list[ab])
                temp_diffncell<-unlist(ncell_difference_filtered_list[ab])
                temp_prevncell<-unlist(prev_ncell_filtered_list[ab])
                temp_newncell<-unlist(new_ncell_filtered_list[ab])#should be aways the same value
                temp_newncell<-temp_newncell[1]
                temp_perc_ncell<-temp_prevncell/temp_newncell
            
                #starting with position of min_dist
                rel_min<-rel_pos_min_dist[ab]
                # changing this minimum eventually
                if ( rel_pos_min_dist[ab] != rel_pos_min_ncell[ab])
                    {
                    dist_mincell<-temp_dist[rel_pos_min_ncell[ab]]
                    dist_mindist<-temp_dist[rel_pos_min_dist[ab]]
                
                    if ( abs(dist_mincell-dist_mindist)<6 & (temp_perc_ncell[rel_pos_min_dist[ab]]>ncell_insteadof_dist_fraction_max | temp_perc_ncell[rel_pos_min_dist[ab]]<ncell_insteadof_dist_fraction_min)) rel_min<-rel_pos_min_ncell[ab]
                    if (abs(dist_mincell-dist_mindist)<1) rel_min<-rel_pos_min_ncell[ab]
                    }
                #find absolute index of chosen rel_min
                temp_absol_pos<-unlist(absol_indexes_filtered_list[ab])
                final_absol_index[ab]<-temp_absol_pos[rel_min]
                if (temp_perc_ncell[rel_min]>ncell_fraction_max_cutout | temp_perc_ncell[rel_min]<ncell_fraction_min_cutout) final_absol_index[ab]<-NA
        
                }
            final_absol_index<-final_absol_index[!is.na(final_absol_index)]
            new_prev_within_dist<-new_prev_within_dist[final_absol_index,]
            new_prev_within_dist<-array(new_prev_within_dist,dim=c(length(new_prev_within_dist)/2,2))
            final_new_prev_withmirror<-new_prev_within_dist
        
            #eliminating duplicated columns
            if (sum(duplicated(final_new_prev_withmirror[,2]))>0)
                {
                cond<-final_new_prev_withmirror[,2] %in% final_new_prev_withmirror[duplicated(final_new_prev_withmirror[,2]),2]
                duplic_prev_eddies<-final_new_prev_withmirror[cond,]
                idf<-factor(duplic_prev_eddies[,2])
                level<-as.integer(levels(idf))
                absol_index<-1:length(duplic_prev_eddies[,2])
            
                dist_duplic_prev_eddies<-sqrt(dist_matr_x[duplic_prev_eddies]^2+dist_matr_y[duplic_prev_eddies]^2)
                rel_min_dist_list<-array(tapply(dist_duplic_prev_eddies,idf,which.min))
                absol_index_list<-array(tapply(absol_index,idf,c))
            
                final_absol_index<-array(0,dim=c(length(level)))
            
                for( ab in 1:length(level))
                    {
                    temporal<-unlist(absol_index_list[ab])
                    final_absol_index[ab]<-temporal[rel_min_dist_list[ab]]
                    }
                eddies_todelete_index<-1:dim(duplic_prev_eddies)[1] #row are not duplicated
            
                eddies_to_delete<-duplic_prev_eddies[!(eddies_todelete_index %in% final_absol_index),]
                eddies_to_delete<-array(eddies_to_delete,dim=c(length(eddies_to_delete)%/%2,2))
                final_new_prev_withmirror<-final_new_prev_withmirror[!(final_new_prev_withmirror[,1] %in% eddies_to_delete[,1]),]
                }
                
            #---final_new_prev_withmirror has got non duplicated eddies, now we need to discern mirror
            dim(final_new_prev_withmirror)<-c(length(final_new_prev_withmirror) %/% 2 ,2)
            final_new_prev<-final_new_prev_withmirror
            found_col_mirror<-final_new_prev_withmirror[final_new_prev[,2]>eddies_tot_num[rel_time-1],]
            found_col_mirror<-array(found_col_mirror,c((length(found_col_mirror) %/% 2),2))
            pos_mirror_in<-which(final_new_prev_withmirror[,2]>eddies_tot_num[rel_time-1])
        
            if (length(pos_mirror_in)==0) {rel_pos_i_b<-NULL}else { rel_pos_i_b<- 1:length(pos_mirror_in)}
            for (i_b in rel_pos_i_b)
                {
                col_mirrored<-mirror_col[ mirror_col[,2] == found_col_mirror[i_b,2],1]
            
                if (sum(final_new_prev_withmirror[,2]==col_mirrored)>0 )
                    {
                    found_col_mirrored<-final_new_prev_withmirror[final_new_prev_withmirror[,2]>eddies_tot_num[rel_time-1],]
                    found_col_mirrored<-array(found_col_mirrored,c((length(found_col_mirrored) %/% 2),2))
                    pos_mirrored_in<-which(final_new_prev_withmirror[,2]==col_mirrored)
                
                    min_mirror<-sqrt(dist_matr_x[found_col_mirror[i_b,1],found_col_mirror[i_b,2]]^2+dist_matr_y[found_col_mirror[i_b,1],found_col_mirror[i_b,2]]^2)
                    min_mirrored<-sqrt(dist_matr_x[found_col_mirrored[1],found_col_mirrored[2]]^2+dist_matr_y[found_col_mirrored[1],found_col_mirrored[2]]^2)
                    if (min_mirror<min_mirrored)
                        {
                        final_new_prev[pos_mirror_in[i_b],2]<-col_mirrored
                        final_new_prev[pos_mirrored_in,]<-c(NA,NA)
                        }
                    else
                        {
                        final_new_prev[pos_mirror_in[i_b],]<-c(NA,NA)
                        }
                    }
                else
                    {
                    final_new_prev[pos_mirror_in[i_b],2]<-col_mirrored
                    }
                }
            final_new_prev<-final_new_prev[!is.na(final_new_prev[,1]),]
            dim(final_new_prev)<-c(length(final_new_prev) %/% 2 ,2)
            index_ed_t[rel_time,final_new_prev[,1]]<-index_ed_t[rel_time-1,final_new_prev[,2]]
        }
            }
    }
    #second run
     if(rel_time>2) index_ed_t<-check_tt2()
     
    write(eddies_tot_num[rel_time], file = file_eddy_index,ncolumns = 1,append = TRUE, sep = " ")
     if(eddies_tot_num[rel_time]>0) write(index_ed_t[rel_time,1:eddies_tot_num[rel_time]], file = file_eddy_index,ncolumns = eddies_tot_num[rel_time],append = TRUE, sep = " ")

    ed_new<- index_ed_t[rel_time,]
    ed_prev<-index_ed_t[rel_time-1,]

    ed_new<-ed_new[!is.na(ed_new)]
    ed_prev<-ed_prev[!is.na(ed_prev)]

    first_track<-ed_new %in% ed_prev
    cat(" tracked first=",sum(first_track))

    if(rel_time>2)
        {
        ed_new<-ed_new[!first_track]
        ed_prev_prev<-index_ed_t[rel_time-2,]
        ed_prev_prev<-ed_prev_prev[!is.na(ed_prev_prev)]
        cat(" sec=",sum(ed_prev_prev %in% ed_new)," ")
        }

    end.time<-Sys.time()
    exec.time<-end.time-start.time
    print(exec.time)
}
