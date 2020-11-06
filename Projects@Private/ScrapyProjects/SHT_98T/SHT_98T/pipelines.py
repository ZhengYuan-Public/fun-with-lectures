# -*- coding: utf-8 -*-

# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://doc.scrapy.org/en/latest/topics/item-pipeline.html
import json


class Sht98TPipeline(object):

    def process_item(self, item, spider):
        contents = json.dumps(dict(item), indent=4, sort_keys=True, ensure_ascii=False)
        with open("./threads_data.json", "ab+") as f:
            f.write(contents.encode("utf-8"))
        return item
