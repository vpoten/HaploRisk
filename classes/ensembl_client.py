import sys
import urllib
import urllib2
import json
import time
import requests


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.max_post = 200

    def perform_rest_action(self, endpoint, hdrs=None, params=None, post_data=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urllib.urlencode(params)

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        if post_data is None:
            return self.__perform_rest_get(endpoint, hdrs, params)

        return self.__perform_rest_post(endpoint, hdrs, post_data)

    def __perform_rest_get(self, endpoint, hdrs, params):
        data = None
        try:
            request = urllib2.Request(self.server + endpoint, headers=hdrs)
            response = urllib2.urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib2.HTTPError, e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    return self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
                return None
        return data

    def __perform_rest_post(self, endpoint, hdrs, post_data):
        r = requests.post(self.server + endpoint, headers=hdrs, data=json.dumps(post_data))
        if not r.ok:
            sys.stderr.write('Request failed for {0}: Status code: {1}\n'.format(endpoint, r.status_code))
            return None
        return r.json()

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None

    def get_snps(self, ids, species='human'):
        ext = "/variation/" + species
        snps = {}
        start = 0

        while start < len(ids):
            result = self.perform_rest_action(ext, hdrs={"Accept": "application/json"},
                                              post_data={'ids': ids[start:start + self.max_post]})
            start += len(result)
            snps.update(result)

        return snps

########################################
# Examples of use
########################################

# def run(species, symbol):
#     client = EnsemblRestClient()
#     variants = client.get_variants(species, symbol)
#     if variants:
#         for v in variants:
#             print '{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v)
#
#
# if __name__ == '__main__':
#     if len(sys.argv) == 3:
#         species, symbol = sys.argv[1:]
#     else:
#         species, symbol = 'human', 'BRAF'
#
#     run(species, symbol)
