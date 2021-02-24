# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# https://doc.scrapy.org/en/latest/topics/items.html

import scrapy


class <removed_name>(scrapy.Item):
    # define the fields for your item here like:
    movie_number_title = scrapy.Field()
    movie_pics_links = scrapy.Field()
    magnet_link = scrapy.Field()
    torrent_link = scrapy.Field()
    torrent_name = scrapy.Field()


class TorrentItem(scrapy.Item):
    file_urls = scrapy.Field()
    files = scrapy.Field()


class ImageItem(scrapy.Item):
    image_urls = scrapy.Field()
    images = scrapy.Field()
