# -*- coding: utf-8 -*-
import scrapy
from ..items import Sht98TItem


class <removed_name>(scrapy.Spider):
    name = 'sht'
    # allowed_domains = ['98rewer.com']
    base_url = '<removed_url>'
    start_urls = [
        '<removed_url>',
                  ]
    page_links = []
    thread_links = []

    def parse(self, response):
        last_page_num = response.xpath("//div[@class='pg']/a[@class='last']/@href").get()
        total_num_pages = int(last_page_num.split('.')[0].split('-')[-1])

        for i in range(total_num_pages):
            page_link = last_page_num.split('.')[0].rsplit('-', 1)[0] + '-' + str(i+1) + '.html'
            page_link = response.urljoin(page_link)
            self.page_links.append(page_link)

        for page in self.page_links:
            yield response.follow(url=page, callback=self.get_threads_in_page)

    def get_threads_in_page(self, response):
        threads_in_page = response.xpath(
            "//tbody[contains(@id,'normalthread')]//td[@class='icn']//a/@href").getall()
        self.thread_links.extend(threads_in_page)

        for thread in self.thread_links:
            # print(f"Parsing thread {thread}")
            yield response.follow(url=thread, callback=self.parse_thread)

    def parse_thread(self, response):
        json_item = Sht98TItem()
        try:
            json_item['movie_number_title'] = response.xpath("//span[@id='thread_subject']/text()").getall()
        except Exception:
            pass
        try:
            json_item['movie_pics_links'] = response.xpath("//td[@class='t_f']//img/@file").getall()
        except Exception:
            pass
        try:
            json_item['magnet_link'] = response.xpath("//div[@class='blockcode']/div//li/text()").getall()
        except Exception:
            pass
        try:
            json_item['torrent_name'] = response.xpath("//p[@class='attnm']/a/text()").getall()
        except Exception:
            pass
        try:
            json_item['torrent_link'] = self.base_url + response.xpath("//p[@class='attnm']/a/@href").getall()[0]
        except Exception:
            pass
        yield json_item

        # try:
        #     torrent_link = self.base_url + response.xpath("//p[@class='attnm']/a/@href").getall()
        #     yield {'file_urls': torrent_link}
        # except Exception:
        #     pass
        # try:
        #     movie_pics_links = response.xpath("//td[@class='t_f']//img/@file").getall()
        #     yield {'image_urls': movie_pics_links}
        # except Exception:
        #     pass




