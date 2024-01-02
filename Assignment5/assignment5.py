#!/usr/bin/env python3

"""python3 assignment1.py -n <aantal_cpus> [OPTIONEEL: -o <output csv file>] fastabestand1.fastq
[fastabestand2.fastq ... fastabestandN.fastq]"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import csv
from pyspark.sql import SparkSession
from pyspark.sql.functions import avg, col
from pyspark.sql.types import DecimalType


def create_session():
    """Create Spark session"""
    spark = (
        SparkSession.builder.master("local[16]")
        .config("spark.executor.memory", "4g")
        .config("spark.driver.memory", "4g")
        .getOrCreate()
    )
    return spark


def read(spark):
    """Run function"""
    df = spark.read.csv(
        "/data/datasets/EBI/interpro/refseq_scan/bacteria.nonredundant_protein.1.protein.faa.tsv",
        sep=r"\t",
        header=False,
    )
    return df


def get_unique_ids(df):
    """Get unique IDS"""
    ids = df.select("_c11").distinct()
    return ids, ids.count()


def get_anno_avg(df):
    """Get anno avg"""
    return df.agg(avg("_c8")).first()[0]


def get_most_occuring(df):
    """Get most occuring"""
    return df.groupBy("_c13").count().sort("count", ascending=False).first()[0]


def get_aa_length(df):
    """Get aa length"""
    return df.select("_c2").first()[0]


def get_top_10(df):
    """Get top 10"""
    return df.groupBy("_c11").count().sort("count", ascending=False).take(10)


def get_top_10_percentage(df):
    """Get top 10 percentage."""
    with_span = df.withColumn("span", (df._c7 - df._c6) / df._c2)
    subset = (
        with_span.sort("span", ascending=False)
        .select(col("_c11"), col("span"), col("_c12"))
        .where(col("span") >= 0.9)
    )
    return subset


def get_top_10_text(df):
    """Get top 10 text."""
    return df.groupBy("_c12").count().sort("count", ascending=False).take(10)


def get_bottom_10_text(df):
    """Get bottom 10."""
    return df.groupBy("_c12").count().sort("count", ascending=True).take(10)


def get_top_10_percentage_text(subset):
    """Get top 10 percentage text"""
    return subset.groupBy("_c12").count().sort("count", ascending=False).take(10)


def get_coefficient(df):
    """Return coefficient between protein size and number of InterPRO annotations"""
    df = df.withColumn("num", df._c8.cast(DecimalType(10, 10)))
    cor = df.withColumn("num2", df._c2.cast(DecimalType(10, 10))).corr("num", "num2")
    return cor


def main():
    """Main function"""
    ses = create_session()
    df = read(ses)
    id_df, id_count = get_unique_ids(df)
    avg = get_anno_avg(df)
    most_occuring = get_most_occuring(df)
    aa_length = get_aa_length(df)
    top_10 = get_top_10(df)
    top_10_percentage_df = get_top_10_percentage(df)
    top_10_percentage = top_10_percentage_df.take(10)
    top_10_text = get_top_10_text(df)
    bottom_10_text = get_bottom_10_text(df)
    top_10_percentage_text = get_top_10_percentage_text(top_10_percentage_df)
    cor = get_coefficient(df)

    with open("output.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "1",
                id_count,
            ]
        )
        writer.writerow(["2", avg])
        writer.writerow(["3", most_occuring])
        writer.writerow(["4", aa_length])
        writer.writerow(["5", top_10])
        writer.writerow(["6", top_10_percentage])
        writer.writerow(["7", top_10_text])
        writer.writerow(["8", bottom_10_text])
        writer.writerow(["9", top_10_percentage_text])
        writer.writerow(["10", cor])

    df.explain(True)


if __name__ == "__main__":
    main()
