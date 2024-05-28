# Assuming `dds` is your DESeq2 object and you have raw counts in `counts(dds)`
raw_counts <- counts(dds)

# Convert to data frame for ggplot2
raw_counts_df <- as.data.frame(raw_counts)
raw_counts_df <- gather(raw_counts_df, key = "sample", value = "count")

# Density plot of raw counts
ggplot(raw_counts_df, aes(x = log2(count + 1), color = sample)) +
  geom_density() +
  ggtitle("Density Plot of Raw Counts") +
  xlab("log2(count + 1)") +
  ylab("Density") +
  theme_minimal()

# Boxplot of raw counts
ggplot(raw_counts_df, aes(x = sample, y = log2(count + 1))) +
  geom_boxplot() +
  ggtitle("Boxplot of Raw Counts") +
  xlab("Sample") +
  ylab("log2(count + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Assuming `rld` is your DESeq2 object after rlog normalization
rlog_counts <- assay(rld)

# Convert to data frame for ggplot2
rlog_counts_df <- as.data.frame(rlog_counts)
rlog_counts_df <- gather(rlog_counts_df, key = "sample", value = "count")

# Density plot of rlog-normalized counts
ggplot(rlog_counts_df, aes(x = count, color = sample)) +
  geom_density() +
  ggtitle("Density Plot of rlog-Normalized Counts") +
  xlab("rlog(count)") +
  ylab("Density") +
  theme_minimal()

# Boxplot of rlog-normalized counts
ggplot(rlog_counts_df, aes(x = sample, y = count)) +
  geom_boxplot() +
  ggtitle("Boxplot of rlog-Normalized Counts") +
  xlab("Sample") +
  ylab("rlog(count)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

