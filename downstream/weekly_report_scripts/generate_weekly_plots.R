#PARSING COMMAND-LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("
  All 4 arguments are required.
  1st argument - full path to sample_summary_stats.csv file
  2nd argument - full path to pipeline summary report file 
    (e.g. 2021-08-31-13-49_ft_65_mutation_report.csv)
  3rd argument - full path to mutation_statistics.csv file
  4th argument - full path to the directory
    where the plot figures should be stored
  ", call.=FALSE)
} 

#LIBRARIES
library(htmlwidgets)
library(plotly)
library(randomcoloR)
library(ggplot2)
library(tidyverse)
library(scales)
library(lubridate)
library(tsibble)


#PALETTE DEFINITIONS
# compressed_palette <-distinctColorPalette(1000) #TO GENERATE NEW PALETTE 
# lapply(compressed_palette, cat, "\n", file="my_palette.txt", append=TRUE) #TO WRITE NEW PALETTE TO A FILE
compressed_palette <- scan("my_palette.txt", what="", sep="\n") #TO READ PALETTE FROM FILE

#PATHS
path <- args[1]#PATH TO SUMMARY FILE
curr_batch_path <- args[2]#PATH TO INHOUSE PIPELINE REPORT
mut_stat_path <- args[3]#PATH TO MUTATION STATISTICS FILE
output_path <- args[4]#PATH TO OUTPUT DIRECTORY
ggsave_width <- 20#PLOT WIDTH IN CM FOR GGPLOT OBJECTS
ggsave_height <- 20#PLOT HEIGHT IN CM FOR GGPLOT OBJECTS

#READING DATA 
excel_data <- read.csv(path,  encoding = "UTF-8")#SUMMARY DATA
batch_data <- read.csv(curr_batch_path, encoding = "UTF-8")#CURRENT BATCH DATA
batch_mut_data <- read.csv(mut_stat_path, encoding = "UTF-8")#MUTATION STATISTICS DATA
# found_data <- read.csv(found_path, encoding = "UTF-8")#FOUND_SAMPLES DATA

#PRE-PROCESSING SUMMARY DATA
excel_data[, 1:4] <- sapply(excel_data[, 1:4], as.character)#CASTING FIRST 4 COLUMNS TO CHARACTER
excel_data <- excel_data %>% filter(!sampling_date == "" & !sampling_date == "0" & !sampling_date == "NaT" & !sampling_date == "44284")#FILTERING INVALID DATA
excel_data$sampling_date <- as.Date(excel_data$sampling_date, format="%Y-%m-%d", origin="1970-01-01")#CASTING DATE FROM CHAR TO DATE.OBJECT
grouped_data <- aggregate(excel_data, by=list(excel_data$sampling_date, excel_data$lineage), FUN=length)[, 1:3]#COUNT NUMBER OF SAMPLES FOR EACH DATE BY LINEAGES
# grouped_data <- data.frame(matrix(unlist(grouped_data), nrow=length(grouped_data), byrow=TRUE))

# grouped <- grouped_data %>% group_by(Group.1) %>% mutate(SamplePer=receiving_lab_sample_id/sum(receiving_lab_sample_id)) %>% ungroup()
# print(grouped_data)
rel_grouped_data <- grouped_data %>% group_by(Group.1) %>% mutate(SamplePer=receiving_lab_sample_id/sum(receiving_lab_sample_id)) %>% ungroup()#COUNT FRACTION OF EACH LINEAGE ON EVERY DATE
rel_grouped_data <- rel_grouped_data %>% mutate(week = yearweek(Group.1))#CONVERT EACH DATE TO WEEK FORMAT - ADD TO NEW COLUMN
rel_grouped_data_week <- rel_grouped_data %>% group_by(week, Group.2) %>% summarise(sum_precip = sum(receiving_lab_sample_id))#SUMMARIZE NUMBER OF SAMPLES PER LINEAGE PER WEEK
rel_grouped_data_week <- as.data.frame(rel_grouped_data_week)#CONVERT LIST TO DATAFRAME
rel_grouped_data_week <- rel_grouped_data_week %>% filter(sum_precip > 0)#KEEP DATES BASED ON NUMBER OF SAMPLES
rel_grouped_data_week_rel <- rel_grouped_data_week %>% group_by(week) %>% mutate(percentage=sum_precip/sum(sum_precip)) %>% ungroup()#COMPUTE FRACTION OF SAMPLES PER WEEK PER LINEAGE
rel_grouped_data_week_rel$week <- as.Date(rel_grouped_data_week_rel$week, "%Y-%U-%u")#CONVERT WEEK COLUMN TO A GIVEN STRING FORMAT FOR PLOTTING
weeks <- unique(rel_grouped_data_week_rel$week)#GET UNIQUE WEEKS
lineages <- unique(rel_grouped_data_week_rel$Group.2)#GET UNIQUE LINEAGES
combinations <- expand.grid(week = weeks, lineages = lineages)#GET ALL COMBINATIONS OF WEEKS && LINEAGES
data <- full_join(rel_grouped_data_week_rel, combinations, by = c("week" = "week", "Group.2" = "lineages")) %>% 
	mutate(percentage = ifelse(is.na(percentage), 0, percentage)) %>% arrange(week, Group.2)#ADD 0 VALUES FOR LINEAGES THAT ARE NOT REPRESENTED ON GIVEN WEEK
data <- data %>% filter(week > as.Date("2021-01-01"))

# #PRE-PROCESSING BATCH DATA
batch_data$receiving_lab_sample_id <- as.character(batch_data$receiving_lab_sample_id)#CASTING SAMPLE ID AS CHARACTER
batch_data$sampling_date <- as.Date(batch_data$sampling_date, format="%Y-%m-%d", origin="1970-01-01")#CASTING DATE STRING AS DATE OBJECT
batch_data$seq_date <- as.Date(batch_data$seq_date, format="%Y-%m-%d", origin="1970-01-01")#CASTING SEQUENCING DATE STRING AS DATE OBJECT
batch_data$covered_reference <- with(batch_data, ceiling(genome_length*(1-genome_N_percentage/100)))#COMPUTING NUMBER OF REFERENCE BASES COVERED X15 OR MORE
batch_data$AVERAGE_COVERAGE <- round(batch_data$AVERAGE_COVERAGE, 2)#ROUNDING-UP AVERAGE COVERAGE TO TWO DIGITS
batch_data$covered_reference_percent <- with(batch_data, round(100*covered_reference*(1/29903),2))#COMPUTING % OF REFERENCE BASES COVERED X15 OR MORE
grp_batch_data <- aggregate(batch_data, by=list(batch_data$sampling_date, batch_data$lineage), FUN=length)[, 1:3]#COUNT NUMBER OF SAMPLES PER SAMPLING DATE PER LINEAGE
# grp_batch_data <- grp_batch_data %>% group_by(Group.1) %>% mutate(SamplePer=receiving_lab_sample_id/sum(receiving_lab_sample_id)) %>% ungroup()#COMPUTE FRACTION OF SAMPLES PER SAMPLING DATE PER LINEAGE
# write.csv(grp_batch_data, paste(output_path, 'test.csv', sep=''))
#PRE-PROCESSING MUTATION DATA
batch_mut_data <- batch_mut_data %>% rename(Percent_of_samples = X._of_samples)#RENAME POORLY NAMED COLUMN
batch_mut_data <- batch_mut_data %>% filter(Percent_of_samples > 3)#FILTER MUTATION DATA - KEEP MUTATIONS THAT OCCUR IN MORE THAT 3 % OF SAMPLE IN THE BATCH
batch_mut_data$Mutation <- factor(batch_mut_data$Mutation, levels = batch_mut_data[order(batch_mut_data$Gene), "Mutation"]) #REORDER MUTATIONS BASED ON GENE
batch_mut_data$Gene <- fct_rev(batch_mut_data$Gene) #REVERSE THE ORDER OF GENES TO MATCH THE ORDER OF MUTATIONS
# batch_mut_data <- batch_mut_data %>% group_by(Gene) %>% mutate(count_gene_occurr = n()) #COUNT NUMBER OF TIMES EACH GENE IS MUTATED
# batch_mut_data$Gene <- reorder(batch_mut_data$Gene, -batch_mut_data$count_gene_occurr) #REORDER BATCH MUTATION DATA IN ORDER FROM LEAST TO MOST MUTATED GENE

#PLOTTING SUMMARY DATA
nedēļa <- data$week
paraugu_daļa <- data$percentage
s_a_plot <- data %>% ggplot( #INITIALIZE PLOT VARIABLE & SUPPLY DATA TO GGPLOT
    aes(
        x = nedēļa, #SUPPLY X AXES DATA
        y = paraugu_daļa)) + #SUPPLY Y AXES DATA
    geom_area( #INITIALIZE GEOMETRIC AREA PLOT
        position='stack', #STACKED LAYOUT
        stat = "identity", #DO NOT SUMMARIZE DATA PRIOR TO PLOTTING, PLOT AS IS
        aes(fill = Group.2), #LEGEND & COLOURING - BY LINEAGES
        color='black', alpha = 0.45) +#AREA OUTLINE COLOUR & OPACITY
    labs(
        x = "", #BLANK X-AXIS LABEL
        y = "Paraugu daļa", #Y-LABEL
        fill = 'Paveidi') + #LEGEND LABEL
    scale_x_date(
        date_breaks = "1 months", #MAYOR BRAKES ON DATE-BASED X-AXIS
        date_minor_breaks = "1 week", #MINOR BRAKES ON DATE-BASED X-AXIS
        date_labels = "%b %Y") + #TICK LABE FORMATE ON DATE-BASED X-AXIS
    theme(
        text = element_text(size=14), #AXIS LABEL & LEGEND FONT SIZE
        axis.text.x = element_text(size=14), #X-AXIS TICK LABEL FONT SIZE
        axis.text.y = element_text(size=14), #Y-AXIS TICK LABEL FONT SIZE
        panel.grid.major.x = element_blank(), #REMOVE VERTICAL GRID LINES
        panel.grid.major.y = element_blank(), #FORMAT HORIZONTAL GRID LINES
        panel.background = element_blank(), #REMOVE BACKGROUND GRID PANEL COLOUR
        axis.line = element_line(colour = "black"),#CHANGE AXIS COLORS
        axis.title.x = element_text(vjust=-1.5),#ADD SPACE BETWEEN TICK LABELS AND X AXIS LABEL
        axis.title.y = element_text(vjust=3),#ADD SPACE BETWEEN TICK LABELS AND Y AXIS LABEL
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + #ADD MARGINS TO THE PLOT
    scale_fill_manual(values=as.vector(compressed_palette)) #
# ggplotly(s_a_plot)
# saveWidget(as_widget(s_a_plot), "1_test.html")
ggsave(
    paste(output_path,'plot7.png', sep=''), #PATH & FILE NAME WHERE THE PLOT WILL BE SAVED
    width = 41.68, #PLOT WIDTH IN CM
    height = 26.73, #PLOT HEIGHT IN CM
    units = "cm") #PLOT DIMENSION UNITS

#PLOTTING ABSOLUTE BATCH DATA
batch_abs_plot <- grp_batch_data %>% ggplot( #INITIALIZE PLOT VARIABLE & SUPPLY DATA TO GGPLOT
    aes(
        y=receiving_lab_sample_id, #SUPPLY Y AXES DATA
        x=Group.1)) + #SUPPLY X AXES DATA
    geom_bar( #INITIALIZE STACKED BAR PLOT
        position="stack", #STACKED LAYOUT
        stat="identity", #PLOT DATA AS IS
        aes(
            fill = forcats::fct_rev(Group.2)), #GENERATE STACKED PLOT WITH BIGGER VALUES ON TOP
            color = 'black', alpha = 0.45) + #OUTLINE COLOUR & OPACITY
        scale_fill_manual(values=as.vector(compressed_palette)) + #ADDING CUSTOM COLOUR PALETTE
    theme(
        axis.text.x = element_text( #FORMAT X-AXIS LABEL TEXT
            angle = 90, #ROTATE COUNTER-CLOCKWISE
            vjust=0.2, #OPTIMAL VERTICAL POSITIONING RELATIVE TO AXIS TICK
            hjust = 0.95), #OPTIMAL HORIZONTAL POSITIONING RELATIVE TO AXIS TICK
        panel.grid.major.x = element_blank(), #REMOVE VERTICAL GRID LINES
        panel.grid.major.y = element_line( #FORMAT VERTICAL GRID LINES
            size=.1, 
            color="black" ),
            panel.background = element_blank(), #REMOVE BACKGROUND GRID PANEL COLOUR
            axis.line = element_line(colour = "black")) + #CHANGE AXIS COLORS
    geom_text(
        aes(label=receiving_lab_sample_id), #ADD PLOTTED VALUE NUMBERS TO STACKED PLOTS
        position=position_stack(reverse=FALSE, vjust=0.5), #SINGLE OPTIMAL LABEL POSITION
        size = 2.2) #LABEL SIZE
plot_4 <- batch_abs_plot + scale_x_date(
    date_breaks = "1 day", #PLOT DATA BY DATE
    date_minor_breaks = "1 day", #PLOT DATA BY DATE
    date_labels = "%e-%m-%Y") + #X-AXIS LABEL FORMAT
    labs(
        x = "Parauga paņemšanas datums", #X-AXIS LABEL
        y = "Paraugu skaits dienā", #Y-AXIS LABEL
        fill = 'Celmi') + #LEGEND LABEL
        theme(
            axis.title.x = element_text(vjust=-1.5), #ADD SPACE BETWEEN TICK LABELS AND X AXIS LABEL
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) #ADD MARGINS TO THE PLOT
ggsave(
    paste(output_path,'plot5.png', sep=''), #PATH & FILE NAME WHERE THE PLOT WILL BE SAVED
    plot=plot_4, #SUPPLYING PLOT OBJECT
    width = 30, #PLOT WIDTH IN CM
    height = 16, #PLOT HEIGHT IN CM
    units = "cm") #PLOT DIMENSION UNITS

#PLOTTING MUTATION DATA
batch_mut_plot <-ggplot( #INITIALIZE PLOT OBJECT
    data=batch_mut_data, #SUPPLYING DATA
    aes(
        y=Number_of_samples, #X-AXIS DATA
        # x=reorder(Mutation, count_gene_occurr))) + #Y-AXIS DATA
        x=Mutation)) +
    geom_bar(
        stat="identity", #PLOT AS-IS
        aes(fill=Gene), #LEGEND BASED ON GENE COLUMN
        width = 0.75, #BAR WIDTH
        color = 'black', alpha = 0.45) + #BAR OUTLINE COLOUR
    scale_fill_manual(values=as.vector(compressed_palette)) +#ADDING CUSTOM COLOUR PALETTE
    theme(
        axis.text.x = element_text(angle = 0, vjust=0.2, hjust = 0.95), #X-AXIS TICK LABEL FORMATTING
        panel.grid.major.y = element_blank(), #REMOVE Y-GRIDLINES
        panel.grid.major.x = element_line(#FORMAT X-GRIDLINES
            size=.1, #GRIDLINE THICKNESS
            color="black" ), #GRIDLINE COLOUR
        panel.background = element_blank(), #REMOVE GRID BACKGROUND
        axis.line = element_line(colour = "black")) #RECOLOR AXES
plot_5<-batch_mut_plot + labs( 
    x = "Mutācijas", #X-AXIS LABEL
    y = "Paraugu skaits", #Y-AXIS LABES
    fill = 'Gēns') +  #LEGEND LABES
    coord_flip() +
    geom_text(size = 2.5, aes(label = Number_of_samples), #DISPLAY VALUE ON TOP OF THE BAR
        vjust = 0.5, #V-ALIGN VALUE
        hjust = -0.2, #H-ALIGN VALUE
        show.legend = FALSE) + #DO NOT INCLUDE LABEL VALUE AS SEPARATE LEGEND ENTRY
        scale_x_discrete(expand = c(0,3)) + #ADD SPACE BETWEEN X AXIS AND PLOT
        theme(
            axis.title.x = element_text(vjust=-1.5),#ADD SPACE BETWEEN TICK LABELS AND X AXIS LABEL
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) #ADD MARGINS TO THE PLOT
ggsave(
    paste(output_path,'plot6.png', sep=''), #PATH & FILE NAME WHERE THE PLOT WILL BE SAVED
    plot=plot_5, #SUPPLYING PLOT OBJECT
    width = ggsave_width, #PLOT WIDTH IN CM
    height = as.integer(nrow(batch_mut_data)/3), #PLOT HEIGHT IN CM
    units = "cm") #PLOT DIMENSION UNITS


#PROCESSING QC METRICS DATA
###FUNCTIONS###
split_df <- function(df) {
    #FUNCTION IS USED TO SPLIT BATCH DATA IN 2 PARTS BASED ON RECEIVING_LAB_SAMPLE_ID
    lower_half <- data.frame(df[df$receiving_lab_sample_id <= median(as.integer(df$receiving_lab_sample_id)),]) #LOWER PART INCLUDES SAMPLES WITH ID LOWER THAT MEDIAN ID
    higher_half <- data.frame(df[df$receiving_lab_sample_id > median(as.integer(df$receiving_lab_sample_id)),]) #HEIGHER PART INCLUDES SAMPLES WITH ID HIGHER THAT MEDIAN ID
    return(list(lower_half,higher_half)) #RETURN LIST OF DATAFRAMES
}

plot_quality_metrics <- function(df, number_list, plot_num, output_path, skip_legend=FALSE){
    #FUNCTION IS USED TO CREATE QC METRIC PLOT BASED ON BATCH_DATA (OR PART OF IT)
    ##METRICS TO PLOT
        total_reads <-df$TOTAL_READS 
        av_coverage <- df$AVERAGE_COVERAGE
        covered_reference <- df$covered_reference_percent
    #PLOT MARGINS
        l_margin <- 15
        r_margin <- 1
        t_margin <- 3
        #DIFFERENT BUTTOM MARGINS FOR 1ST AND LATER PLOTS
            if (skip_legend==TRUE){
                b_margin=6
            }else{
                b_margin=12
            }
    png(file=paste(output_path, 'plot', plot_num,'.png', sep=''),width=1500, height=500) #SAVE PLOT AS PNG
    par(mar=c(b_margin, l_margin, t_margin, r_margin) + 0.1, cex.axis=1.3)#INITIALIZE A LAYOUT WITH DEFINED MARGINS, SCALE AXES FONT SIZE UP 30%
    #PLOT_1
        plot( #INITIALIZE PLOT OBJECT
            c(min(number_list):max(number_list)), #DEFINE X-AXIS RANGE
            las=1, #ORIENT TICK LABELS HORIZONTALLY
            axes=F, #ALLOW MULTIPLE Y AXES IN ONE LAYOUT
            covered_reference, #SUPPLY Y DATA  
            xlim = c( #DEFINE LIMIT VALUES TO DISPLAY ON X-AXIS
                as.integer(min(number_list)),
                as.integer(max(number_list))), 
            ylim=c(0,max(covered_reference)), #DEFINE LIMIT VALUES TO DISPLAY ON Y-AXIS
            xlab="", #X-AXIS LABEL
            ylab="", #Y-AXIS LABEL
            type="l", #TYPE OF THE PLOT TO DRAW - LINE PLOT
            main="", #LEGEND LABEL
            lwd=3, #PLOT LINE WIDTH
            col="darkgreen") #PLOT LINE COLOUR
        segments( #DRAW LINES FROM X-AXIS TO DATA POINTS
            x0=c(min(number_list):max(number_list)), #LINE STARTING LEVEL X-COORDINATE
            y0=0, #LINE STARTING LEVEL Y-COORDINATE
            y1=covered_reference, #LINE ENDPOINT VALUES
            col=8) #LINE TYPE (E.G. DASHED LINES)
        points( #ADD POINTS TO DATA ON PLOT
            number_list, #X-AXIS VALUES
            covered_reference, #Y-AXIS VALUES
            pch=20, #POINT TYPES (E.G. CIRCLE, SQUARE ETC.)
            col="darkgreen", #POINT COLOUR
            cex=2.5) #SCALE POINT SIZE 150% UP
        axis(2, #Y-AXIS LOCATION; 2-LEFT
            ylim=c(min(covered_reference),max(covered_reference)), #Y-AXIS DATA RANGE
            col="darkgreen", #Y-AXIS COLOUR
            lwd=2) #Y-AXIS LINE WIDTH
        mtext(2,
            text="Reference pārklāta (min.x15 %)", #Y-AXIS LABEL
            line=2.5, #PADDING FOR Y-AXIS LABEL
            cex=1.6) #SCALE X-AXIS LABEL FONT SIZE 60% UP

    #PLOT_2
        par(
            new=T, #ADD NEW PLOT TO THE LAYOUT
            cex.axis=1.3) #SCALE AXIS TICK LABEL FONT SIZE UP 30%
        plot(c(min(number_list):max(number_list)), 
            las=1,
            axes=F, 
            av_coverage, 
            ylim=c(0,max(av_coverage)), 
            xlab="", 
            ylab="", 
            xlim = c(
                as.integer(min(number_list)),
                as.integer(max(number_list))), 
            type="l",
            lty=2, #PLOT LINE TYPE
            main="",
            lwd=3, 
            col="blue")
        segments(
            x0=c(min(number_list):max(number_list)), 
            y0=0,
            y1=av_coverage, 
            col=8)
        axis(2, 
            ylim=c(
                min(av_coverage),
                max(av_coverage)),
            lwd=2,
            line=4.2, 
            col="blue")
        points(
            number_list, 
            av_coverage,
            pch=20, 
            col="blue",
            cex=2.5)
        mtext(2,
            text="Vidējais pārklājums",
            line=6.5,
            cex=1.6)

    #PLOT_3
        par(new=T, 
            cex.axis=1.3)
        plot(c(min(number_list):max(number_list)), 
            las=1,
            axes=F, 
            total_reads, 
            ylim=c(0,max(total_reads)), 
            xlab="", 
            ylab="", 
            xlim = c(
                as.integer(min(number_list)),
                as.integer(max(number_list))), 
            type="l",
            lty=3, 
            main="",
            lwd=3, 
            col="red")
        segments(
            x0=c(min(number_list):max(number_list)), 
            y0=0,
            y1=total_reads, 
            col=8)
        axis(2, 
        ylim=c(
            min(total_reads),
            max(total_reads)), 
            lwd=2,
            line=8.2, 
            col="red")
        points(
            number_list, 
            total_reads, 
            pch=20, 
            col="red", 
            cex=2.5)
        mtext(2,
            text="Nolasījumu skaits",
            line=10.5, 
            cex=1.6)
        axis(1, 
            at=min(number_list):max(number_list), 
            labels=number_list)
        if (skip_legend==FALSE){    #TO NOT ADD LEGEND IF IT IS NOT NEEDED
            legend(
                "bottom",  #LEGEND POSITION
                xpd=TRUE, #DRAW LEGEND OUTSIDE OF PLOT AREA
                inset = c(0.2,-0.6), #LEGEND POSITIONING
                cex=1.6, #SCALE LEGEND SIZE UP 60%S
                legend=c( #LEGEND LABELS
                    "Reference pārklāta (min.x15 %)", 
                    "Vidējais pārklājums",
                    "Nolasījumu skaits"),
                lty=c(1,2,3)) #LINE TYPES TO DISPLAY ON LEGEND
            mtext(
                "Parauga ID",
                side=1,
                col="black", 
                line=3.2,
                cex=1.6)
        }else{
            mtext(
                "Parauga ID",
                side=1,
                col="black", 
                line=3.5,
                cex=1.6)
        }
    dev.off() #SAVE PLOT
}

#PRE-PROCESSING QC METRIC
chunk_len <- nrow(batch_data)/4 + 1 #DEVIDE NUMBER OF SAMPLES BY NUMBER OF REQUIRED CHUNKS TO GET OPTIMAL CHUNK SIZE (+1 TO AVOID INCORRECT SPLITTING FOR ODD NUMBER OF SAMPLES)
target <- batch_data$receiving_lab_sample_id #GETTING ORDER OF SAMPLES IN FOUND SAMPLES REPORT
batch_data <- batch_data[match(target, batch_data$receiving_lab_sample_id),] #REORDERING SAMPLES ACCORDING TO FOUND_SAMPLES REPORT
batch_data$receiving_lab_sample_id <- seq.int(nrow(batch_data)) #ADD NUMERIC ID COLUMN TO SPLIT UPON
write_csv(batch_data, 'test.csv')
# write.csv(batch_data, paste(output_path, 'test.csv', sep=''))
chunk_list <- list(batch_data) #INITIALIZE LIST OF CHUNKS WITH STARTING DATAFRAME
while (nrow(data.frame(chunk_list[1])) > chunk_len) { #UNTIL CHUNK SIZE IS NOT SMALLER THAN CHUNK_LEN VALUE
    temp_list <- list() #EMPTY LIST TO STORE CURRENT CHUNKS
    for (chunk in chunk_list) { #ITERATE OVER CURRENT CHUNKS
        split_result = split_df(chunk) #SPLIT EACH CHUNK IN TWO
        temp_list <- append(temp_list, split_result[1]) #ADD CHUNK WITH LOWER IDS FIRST
        temp_list <- append(temp_list, split_result[2]) #ADD CHUNK WITH HIGHER IDS SECOND

    }
    chunk_list <- temp_list #WHEN ALL CHUNKS WERE SPLIT, THE SPLIT PRODUCTS IS STORED IN TEMP_LIST - ASSIGNED TO CHUNK_LIST
}
#REPORTING DATA INTEGRITY STATUS
total_df = data.frame() #TO STORE DATA AS IT IS RESTORED FROM CHUNKS
for (chunk in chunk_list){ #MERGE CHUNKS FROM CHUNK LIST INTO DATAFRAME
    print(min(chunk$receiving_lab_sample_id))#VERIFYING THAT CHUNKS ARE STORED IN ORDER OF INCREASING ID - PRINT LOWEST ID IN CHUNK
    print(max(chunk$receiving_lab_sample_id))#VERIFYING THAT CHUNKS ARE STORED IN ORDER OF INCREASING ID - PRINT HIGHEST ID IN CHUNK
    total_df = rbind(total_df, data.frame(chunk)) #MERGE CHUNK INTO DATAFRAME
}
paste("Split batch data into", length(chunk_list), "chunks.") #PRINT NUMBER OF CHUNKS INTO STANDARD OUTPUT
paste(nrow(total_df),"/",nrow(batch_data)) #COMPARE NUMBER OF ROWS IN ORIGINAL DATAFRAME TO THE NUMBER OF ROWS IN DATAFRAME RESTORED FROM CHUNKS


#PLOTTING QC METRICS
for (i in 1:length(chunk_list)){ #FOR EACH CHUNK
    skip_legend <- FALSE #ONLY FOR THE FIRST PLOT
    if (i==1){
        number_list <- c(1:nrow(data.frame(chunk_list[i]))) #FOR THE FIRST CHUNK SAMPLES ARE NUMBERED FROM 1
        last_id <- nrow(data.frame(chunk_list[i])) #LAST ID TO BE USED WHEN INCREMENTING SAMPLE NUMBER IN MULTIPLE PLOTS
    }else{
        skip_legend <- TRUE #FOR ALL PLOTS EXCEPT THE FIRST
        number_list <- c((last_id+1):(last_id + nrow(data.frame(chunk_list[i])))) #FOR EVERY SUBSEQUENT PLOT SAMPLES ARE NUMBERED FROM LAST_ID + 1
        last_id <- last_id + nrow(data.frame(chunk_list[i])) #UPDATING LAST ID 
    }
    chunk_df <- data.frame(chunk_list[i]) #GETTING CHUNK DATA
    chunk_df <- chunk_df[order(chunk_df$receiving_lab_sample_id),] #ORDERING CHUNK DATA FROM SMALLEST TO GREATEST RECEIVING_LAB_SAMPLE_ID 
    plot_quality_metrics(chunk_df, number_list,i,output_path, skip_legend) #PLOT CHUNK DATA
}


#PRE-PROCESSING DISTRICT DATA
batch_data <- batch_data %>% filter(!district == "nan", !district == "0")
grp_district_data <- aggregate(batch_data, by=list(batch_data$district, batch_data$lineage), FUN=length)[, 1:3] #COUNT NUMBER OF SAMPLES PER DISTRICT PER LINEAGE
# grp_district_data <- grp_district_data %>% filter(!receiving_lab_sample_id < 2) #REMOVE LINEAGES IF COUNT IS LESS THAN % OF THE TOTAL COUNT
grp_district_data <- grp_district_data %>% group_by(Group.1) %>% mutate(SamplePer=receiving_lab_sample_id/sum(receiving_lab_sample_id)) %>% ungroup()#COUNT FRACTION OF EACH LINEAGE FROM EVERY DISTRICT
grp_district_data <- grp_district_data %>% group_by(Group.1) %>% mutate(SumPer=sum(receiving_lab_sample_id)) %>% ungroup() #TOTAL SAMPLES PER CITY

#PLOTTING LINEAGES BY DISTRICT
paraugu_daļa <- round(grp_district_data$SamplePer,3)
novads <- fct_reorder(grp_district_data$Group.1, desc(grp_district_data$SumPer))
paveids <- forcats::fct_rev(grp_district_data$Group.2)
paraugu_kopā <- grp_district_data$SumPer
paraugu_skaits <- grp_district_data$receiving_lab_sample_id


district_plot <- grp_district_data %>% ggplot( #INITIALIZE PLOT VARIABLE & SUPPLY DATA TO GGPLOT
    aes(
        y=paraugu_daļa, #SUPPLY Y AXES DATA
        x=novads)) + #SUPPLY X AXES DATA (SORTED BY DECR.ORDER OF TOTAL NUMBER OF SAMPLES IN EACH DISTRICT)
    geom_bar( #INITIALIZE STACKED BAR PLOT
        position="stack", #STACKED LAYOUT
        stat="identity", #PLOT DATA AS IS
        aes(
            fill = paveids), #GENERATE STACKED PLOT WITH BIGGER VALUES ON TOP
            color = 'black', alpha = 0.45) + #SUB-BAR OUTLINE COLOUR + OPACITY REDUCED
        scale_fill_manual(values=as.vector(compressed_palette)) +#ADDING CUSTOM COLOUR PALETTE
    theme(
        axis.text.x = element_text( #FORMAT X-AXIS LABEL TEXT
            angle = 90, #ROTATE COUNTER-CLOCKWISE
            vjust=0.2, #OPTIMAL VERTICAL POSITIONING RELATIVE TO AXIS TICK
            hjust = 0.95), #OPTIMAL HORIZONTAL POSITIONING RELATIVE TO AXIS TICK
        panel.grid.major.x = element_blank(), #REMOVE VERTICAL GRID LINES
        panel.grid.major.y = element_line( #FORMAT VERTICAL GRID LINES
            size=.1, 
            color="black" ),
            panel.background = element_blank(), #REMOVE BACKGROUND GRID PANEL COLOUR
            axis.line = element_line(colour = "black")) + #CHANGE AXIS COLORS
    geom_text(
        aes(label=paraugu_skaits), #ADD PLOTTED VALUE NUMBERS TO STACKED PLOTS
        position=position_stack(reverse=FALSE, vjust=0.5), #SINGLE OPTIMAL LABEL POSITION
        size = 2.2) #LABEL SIZE
plot_8 <- district_plot + 
    labs(
        x = "Parauga paņemšanas novads", #X-AXIS LABEL
        y = "Paraugu skaits",
        fill = 'Celmi') + #LEGEND LABEL
        theme(
            axis.title.x = element_text(vjust=0.5), #ADD SPACE BETWEEN TICK LABELS AND X AXIS LABEL
            axis.title.y = element_text(vjust=1.5), #ADD SPACE BETWEEN TICK LABELS AND Y AXIS LABEL
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +#ADD MARGINS TO THE PLOT
        geom_text(
            size = 2.5, 
            data=grp_district_data,
            aes(x = novads, y=1, label=paraugu_kopā),
            position = position_dodge(0.9),
            vjust=0.25, 
            hjust=-4,
            angle = 90)
# plot_8 <- ggplotly(plot_8)
# saveWidget(as_widget(plot_8), "8_test.html")



# ggsave(
#     paste(output_path,'plot8.png', sep=''), #PATH & FILE NAME WHERE THE PLOT WILL BE SAVED
#     plot=plot_8, #SUPPLYING PLOT OBJECT
#     width = ggsave_width, #PLOT WIDTH IN CM
#     height = ggsave_height, #PLOT HEIGHT IN CM
#     units = "cm") #PLOT DIMENSION UNITS