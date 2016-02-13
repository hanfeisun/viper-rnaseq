'use strict';

/* App Module */

var viperApp = angular.module('viperApp', [
  'ui.router',
  'viperAppControllers'
]);


viperApp.config(function($stateProvider, $urlRouterProvider) {

    $urlRouterProvider.otherwise('/');

    $stateProvider

        .state('/', {
            url: '/',
            templateUrl: 'app/partials/full-report.html'
        })
        .state('/home', {
            url: '/home',
            templateUrl: 'app/partials/test.html'
        });
});


