(function() {
  var slice = [].slice;

  define(["underscore", "./console", "./dom", "./events"], function(_, console, dom, events) {
    var addStylesheets, exports, invokeInitializer, isIE, isOpera, loadLibrariesAndInitialize, pathPrefix, rebuildURL, rebuildURLOnIE;
    pathPrefix = null;
    isOpera = Object.prototype.toString.call(window.opera) === '[object Opera]';
    isIE = !!window.attachEvent && !isOpera;
    rebuildURL = function(path) {
      var l;
      if (path.match(/^https?:/)) {
        return path;
      }
      if (!pathPrefix) {
        l = window.location;
        pathPrefix = l.protocol + "//" + l.host;
      }
      return pathPrefix + path;
    };
    rebuildURLOnIE = isIE ? rebuildURL : _.identity;
    addStylesheets = function(newStylesheets) {
      var head, insertionPoint, loaded;
      if (!newStylesheets) {
        return;
      }
      loaded = _.chain(document.styleSheets).pluck("href").without("").without(null).map(rebuildURLOnIE);
      insertionPoint = _.find(document.styleSheets, function(ss) {
        var parent;
        parent = ss.ownerNode || ss.owningElement;
        return parent.rel === "stylesheet ajax-insertion-point";
      });
      head = document.head || document.getElementsByTagName("head")[0];
      _.chain(newStylesheets).map(function(ss) {
        return {
          href: rebuildURL(ss.href),
          media: ss.media
        };
      }).reject(function(ss) {
        return loaded.contains(ss.href).value();
      }).each(function(ss) {
        var element;
        element = document.createElement("link");
        element.setAttribute("type", "text/css");
        element.setAttribute("rel", "stylesheet");
        element.setAttribute("href", ss.href);
        if (ss.media) {
          element.setAttribute("media", ss.media);
        }
        if (insertionPoint) {
          head.insertBefore(element, insertionPoint.ownerNode);
        } else {
          head.appendChild(element);
        }
        return console.debug("Added stylesheet " + ss.href);
      });
    };
    invokeInitializer = function(tracker, qualifiedName, initArguments) {
      var functionName, moduleName, ref;
      ref = qualifiedName.split(':'), moduleName = ref[0], functionName = ref[1];
      return require([moduleName], function(moduleLib) {
        var arg, argsString, fn;
        try {
          if (!functionName && initArguments.length === 0 && !_.isFunction(moduleLib)) {
            console.debug("Loaded module " + moduleName);
            return;
          }
          fn = functionName != null ? moduleLib[functionName] : moduleLib;
          if (fn == null) {
            throw new Error("Could not locate function `" + qualifiedName + "'.");
          }
          if (console.debugEnabled) {
            argsString = ((function() {
              var i, len, results;
              results = [];
              for (i = 0, len = initArguments.length; i < len; i++) {
                arg = initArguments[i];
                results.push(JSON.stringify(arg));
              }
              return results;
            })()).join(", ");
            console.debug("Invoking " + qualifiedName + "(" + argsString + ")");
          }
          return fn.apply(null, initArguments);
        } finally {
          tracker();
        }
      });
    };
    loadLibrariesAndInitialize = function(libraries, inits) {
      console.debug("Loading " + ((libraries != null ? libraries.length : void 0) || 0) + " libraries");
      return exports.loadLibraries(libraries, function() {
        return exports.initialize(inits, function() {
          var i, len, mask, ref, results;
          dom.body.attr("data-page-initialized", "true");
          ref = dom.body.find(".pageloading-mask");
          results = [];
          for (i = 0, len = ref.length; i < len; i++) {
            mask = ref[i];
            results.push(mask.remove());
          }
          return results;
        });
      });
    };
    return exports = _.extend(loadLibrariesAndInitialize, {
      initialize: function(inits, callback) {
        var callbackCountdown, i, init, initArguments, len, qualifiedName, tracker;
        if (inits == null) {
          inits = [];
        }
        console.debug("Executing " + inits.length + " inits");
        callbackCountdown = inits.length + 1;
        tracker = function() {
          callbackCountdown--;
          if (callbackCountdown === 0) {
            console.debug("All inits executed");
            if (callback) {
              return callback();
            }
          }
        };
        for (i = 0, len = inits.length; i < len; i++) {
          init = inits[i];
          if (_.isString(init)) {
            invokeInitializer(tracker, init, []);
          } else {
            qualifiedName = init[0], initArguments = 2 <= init.length ? slice.call(init, 1) : [];
            invokeInitializer(tracker, qualifiedName, initArguments);
          }
        }
        return tracker();
      },
      loadLibraries: function(libraries, callback) {
        var finalCallback, reducer;
        reducer = function(callback, library) {
          return function() {
            console.debug("Loading library " + library);
            return require([library], callback);
          };
        };
        finalCallback = _.reduceRight(libraries, reducer, callback);
        return finalCallback.call(null);
      },
      evalJavaScript: function(js) {
        console.debug("Evaluating: " + js);
        return eval(js);
      },
      focus: function(fieldId) {
        var field;
        field = dom(fieldId);
        if (field) {
          return _.delay((function() {
            return field.focus();
          }), 125);
        }
      },
      handlePartialPageRenderResponse: function(response, callback) {
        var partial, responseJSON;
        responseJSON = response.json || {};
        partial = responseJSON._tapestry;
        delete responseJSON._tapestry;
        if (partial != null ? partial.redirectURL : void 0) {
          if (window.location.href === partial.redirectURL) {
            window.location.reload(true);
          } else {
            window.location.href = partial.redirectURL;
          }
          return;
        }
        addStylesheets(partial != null ? partial.stylesheets : void 0);
        exports.loadLibraries(partial != null ? partial.libraries : void 0, function() {
          _(partial != null ? partial.content : void 0).each(function(arg1) {
            var content, id, zone;
            id = arg1[0], content = arg1[1];
            console.debug("Updating content for zone " + id);
            zone = dom.wrap(id);
            if (zone) {
              return zone.trigger(events.zone.update, {
                content: content
              });
            }
          });
          callback && callback(response);
          return exports.initialize(partial != null ? partial.inits : void 0);
        });
      }
    });
  });

}).call(this);
